#include "compatibility/SF_compatibility.h"
#include "infrastructure/io/mesh/SF_sfmMesh.h"
#include <filesystem>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <set>
#include <cmath>
namespace SF {
namespace {
using P=Model::Parameters;

P expandBoundaryDefault(P value,const std::string& caseDir,const std::vector<std::string>& meshFiles) {
    if(!value.is_object()||!value.contains("default"))return value;
    const auto fallback=value.at("default");value.erase("default");
    std::set<std::string> names;
    SFMMeshIO mesh(caseDir);
    for(const auto& file:meshFiles) {
        SFMMeshBlock block;
        if(mesh.readSingleBlock(file,block))for(const auto& [name,points]:block.pointSets){ (void)points; names.insert(name); }
        std::map<std::string,std::vector<int>> sets;
        if(mesh.readSetsOnly(file,sets))for(const auto& [name,points]:sets){ (void)points; names.insert(name); }
    }
    if(names.empty())throw std::runtime_error("boundary default requires readable mesh patch sets");
    for(const auto& name:names)if(!value.contains(name))value[name]=fallback;
    return value;
}
// Native coefficients are mandatory; only historical input normalization uses defaults.
ThermophysicalContract parseThermophysical(const Model::ObjectDescriptor& object) {
    const auto& input=object.parameters;
    const std::string source=object.source.empty()?object.category+"/"+object.name:object.source;
    Model::Schema{{{"thermoDynamics","object",true},{"properties","object",true}},false}.validate(input,source);
    const auto& selection=input.at("thermoDynamics");
    if (selection.contains("energy")) throw std::runtime_error(source+": energy is WHAT; select it in equations registry.");
    Model::Schema{{{"equationOfState","string",true},{"thermo","string"},{"transport","string"}},false}.validate(selection,source+"/thermoDynamics");
    ThermophysicalContract t;t.native=true;t.source=source;
    t.selection.equationOfState=selection.at("equationOfState").get<std::string>();
    t.selection.thermo=selection.value("thermo",std::string());
    t.selection.transport=selection.value("transport",std::string());
    const auto& props=input.at("properties");
    Model::Schema{{{"equationOfState","object",true},{"thermo","object"},{"transport","object"}},false}.validate(props,source+"/properties");
    const auto number=[&](const P& p,const std::string& name,const std::string& group) {
        if(!p.contains(name)||!p.at(name).is_number())throw std::runtime_error(source+": missing required numeric properties."+group+"."+name);
        const double v=p.at(name).get<double>();
        if(!std::isfinite(v))throw std::runtime_error(source+": nonfinite properties."+group+"."+name);
        t.provenance[name]=source+"/properties."+group+"."+name;
        return v;
    };
    const auto& eos=props.at("equationOfState");
    if(t.selection.equationOfState=="perfectGas") {
        Model::Schema{{{"gamma","number",true},{"R","number",true}},false}.validate(eos,source+"/properties.equationOfState");
        t.gamma=number(eos,"gamma","equationOfState");t.gasConstant=number(eos,"R","equationOfState");
    } else if(t.selection.equationOfState=="rhoConst") {
        Model::Schema{{{"rho","number",true}},false}.validate(eos,source+"/properties.equationOfState");
        t.selection.constantDensity=number(eos,"rho","equationOfState");
    } else throw std::runtime_error(source+": unknown equationOfState '"+t.selection.equationOfState+"'.");
    if(!t.selection.transport.empty()) {
        if(!props.contains("transport"))throw std::runtime_error(source+": missing required properties.transport (mu/Pr).");
        const auto& tr=props.at("transport");
        Model::Schema{{{"mu","number",true},{"Pr","number",true}},false}.validate(tr,source+"/properties.transport");
        t.dynamicViscosity=number(tr,"mu","transport");t.prandtl=number(tr,"Pr","transport");
    } else if(props.contains("transport"))throw std::runtime_error(source+": properties.transport has no selected transport model.");
    t.validate();
    if(props.contains("thermo")) {
        const auto& caloric=props.at("thermo");
        Model::Schema{{{"cp","number"},{"cv","number"}},false}.validate(caloric,source+"/properties.thermo");
        if(t.selection.equationOfState!="perfectGas"||t.selection.thermo!="hConst")
            throw std::runtime_error(source+": explicit cp/cv validation requires perfectGas+hConst.");
        for(const auto* name:{"cp","cv"})if(caloric.contains(name)) {
            const double v=number(caloric,name,"thermo");const double expected=std::string(name)=="cp"?t.cp():t.cv();
            if(v<=0.||std::abs(v-expected)>1.e-12*std::max(1.,std::abs(expected)))
                throw std::runtime_error(source+": inconsistent cp/cv with gamma/R; cp-cv=R and gamma=cp/cv are required.");
        }
    }
    return t;
}

// Typed native sections are stored in CaseAdapter::sections_;
// native decoding lives in SF_nativeDecode.cpp.
}
CaseConfig CaseAdapter::build(const Model::Description& m) {
    native_=true;sections_={};
    EquationCompositionConfig composition;
    std::optional<ThermophysicalContract> thermophysical;
    P nativeAlgorithmParameters=P::object();
    std::optional<System::ExecutionProgram> authoredExecution;
    std::vector<System::NumericalBinding> authoredNumerics;
    bool numericalBindingsDeclared=false;
    // 语义角色 → parameters。同一个角色只能被一个对象占用。
    std::map<std::string,P> roles;
    auto assignRole=[&](const std::string& role,const P& parameters) {
        if(!roles.emplace(role,parameters).second)
            throw std::runtime_error("Multiple model objects target the same built-in configuration: "+role);
    };
    Model::Schema{{{"startTime","number"},{"endTime","number"},{"endStep","integer"},
        {"CFL","number"},{"maxDeltaT","number"},{"writeControl","string",true},
        {"writeInterval","number",true},{"writeInitial","boolean"},{"createMesh","boolean"},
        {"parallel","object"},{"output","object",true}},false}.validate(m.runtime,"runtime.yaml");
    // 每个语义对象直接落到自己的 typed slot；不经过任何磁盘布局名字。
    sections_.runtime=m.runtime;
    sections_.output=m.runtime.at("output");sections_.hasOutput=true;
    if(m.runtime.contains("parallel")) { sections_.parallel=m.runtime.at("parallel");sections_.hasParallel=true; }
    sections_.algorithm=m.solver;sections_.hasAlgorithm=!m.solver.is_null();
    const bool hasTermRecipes=m.numerics.contains("terms");
    if(hasTermRecipes&&(m.numerics.contains("convection")||m.numerics.contains("diffusion")))
        throw std::runtime_error("numerics.yaml cannot mix terms with legacy convection/diffusion selection");
    for(auto it=m.numerics.begin();it!=m.numerics.end();++it) {
        const bool known=it.key()=="time"||it.key()=="convection"
            ||it.key()=="diffusion"||it.key()=="sources"
            ||it.key()=="transport"||it.key()=="pressureCorrection"
            ||it.key()=="terms";
        if(!known)throw std::runtime_error("numerics.yaml: unknown operator group "+it.key());
        if(it.key()=="terms") {
            if(!it.value().is_object())throw std::runtime_error("numerics.yaml: terms must be a mapping");
            for(auto term=it.value().begin();term!=it.value().end();++term)
                if(term.key()!="convection"&&term.key()!="diffusion")
                    throw std::runtime_error("numerics.yaml: unknown term recipe role "+term.key());
        }
    }
    sections_.numerics=m.numerics;
    Model::Schema{{{"files","array",true},{"generator","string"}},false}.validate(m.mesh,"mesh/mesh.yaml");
    meshFiles_=m.mesh.at("files").get<std::vector<std::string>>();
    if(m.runtime.value("createMesh",false)) {
        // Geometry format adapter: blockMeshDict is a mesh format, not the model case format.
        const std::string generator=m.mesh.at("generator");
        std::ifstream input((std::filesystem::path(caseDir_)/generator).string());
        if(!input)throw std::runtime_error("Missing mesh generator: "+generator);
        meshGenerator_=generator;
    }
    using Context=std::function<void(const std::string&,const P&)>;
    Model::FactoryRegistry<Context> factories;
    factories.add("stateRegistry",{{{"use","array",true}},false},
        [&](const Model::ObjectDescriptor& object,Context&) {
            if (composition.stateDeclared) throw std::runtime_error("Multiple solution STATE selections.");
            composition.declared=true;
            composition.stateDeclared=true;
            composition.solutionVariables=object.parameters.at("use").get<std::vector<std::string>>();
            std::set<std::string> selected;
            for (const auto& id:composition.solutionVariables) {
                if (id.empty()) throw std::runtime_error("Empty solution STATE symbol.");
                if (!selected.insert(id).second) throw std::runtime_error("Duplicate solution STATE: "+id);
            }
            composition.stateSelectionOrigin="native stateRegistry";
        });
    factories.add("equationRegistry",{{{"use","array",true},{"add","object"},
        {"extend","object"},{"replace","object"},{"disable","array"}},false},
        [&](const Model::ObjectDescriptor& object,Context&) {
            composition.declared = true;
            if (!composition.equations.empty()) throw std::runtime_error("Multiple WHAT equation selections.");
            composition.equations = object.parameters.at("use")
                .get<std::vector<std::string>>();
            for (const auto& equation : composition.equations) {
                if (equation == "Momentum" || equation == "Continuity"
                    || equation == "Energy" || equation=="PressureConstraint") continue;
                if (equation.find('/') != std::string::npos) {
                    throw std::runtime_error(
                        "Custom equation '"+equation
                        +"' is registered but its Equation compiler is not implemented.");
                }
                throw std::runtime_error("Unknown built-in equation '"+equation+"'.");
            }
            if (object.parameters.contains("add") || object.parameters.contains("extend")
                || object.parameters.contains("replace") || object.parameters.contains("disable")) {
                throw std::runtime_error(
                    "Equation modifications are declared by the input contract, "
                    "but lowering custom modifications is not implemented.");
            }
        });
    factories.add("executionProgram",{{{"root","object",true}},false},
        [&](const Model::ObjectDescriptor& object,Context&) {
            if (authoredExecution) throw std::runtime_error("Multiple authored HOW programs.");
            const auto parse=[&](const auto& self,const P& value)->System::ExecutionScope {
                Model::Schema{{{"kind","string",true},{"id","string"},{"order","integer"},
                    {"equation","string"},{"target","string"},{"targetKind","string"},{"occurrence","string"},
                    {"children","array"},{"repetitions","integer"},{"minimumIterations","integer"},
                    {"terminationSignal","string"}},false}.validate(value,"HOW node");
                System::ExecutionScope node;
                const auto kind=value.at("kind").get<std::string>();
                if (kind=="Sequence") node.kind=System::ExecutionKind::Sequence;
                else if (kind=="Loop") node.kind=System::ExecutionKind::Loop;
                else if (kind=="StageLoop") node.kind=System::ExecutionKind::StageLoop;
                else if (kind=="Commit") node.kind=System::ExecutionKind::Commit;
                else if (kind=="EquationCall") node.kind=System::ExecutionKind::EquationCall;
                else throw std::runtime_error("Unknown HOW node kind: "+kind);
                node.id=value.value("id",std::string());node.order=value.value("order",0);
                node.origin={System::OriginKind::User,"native executionProgram"};
                node.repetitions=value.value("repetitions",1);node.minimumIterations=value.value("minimumIterations",1);
                node.terminationSignal=value.value("terminationSignal",std::string());
                if (node.kind==System::ExecutionKind::EquationCall) {
                    if (!value.contains("equation") || !value.contains("target") || value.contains("children"))
                        throw std::runtime_error("HOW EquationCall requires equation/target and no child scope.");
                    node.step={value.at("equation").get<std::string>(),
                        System::targetFromSyntax(value.at("target").get<std::string>()),value.value("occurrence",std::string())};
                    if (value.contains("targetKind")) {
                        if (node.step.target.kind!=System::TargetKind::Physical)
                            throw std::runtime_error("HOW targetKind cannot override a qualified target.");
                        const auto targetKind=value.at("targetKind").get<std::string>();
                        if (targetKind=="Physical") node.step.target.kind=System::TargetKind::Physical;
                        else if (targetKind=="Working") node.step.target.kind=System::TargetKind::Working;
                        else if (targetKind=="Correction") node.step.target.kind=System::TargetKind::Correction;
                        else if (targetKind=="Workspace") node.step.target.kind=System::TargetKind::Workspace;
                        else throw std::runtime_error("Unknown HOW target kind: "+targetKind);
                    }
                } else {
                    if (value.contains("equation") || value.contains("target") || value.contains("targetKind") || value.contains("occurrence"))
                        throw std::runtime_error("HOW control scope cannot carry an equation target.");
                    if (value.contains("children")) for (const auto& child:value.at("children")) node.children.push_back(self(self,child));
                }
                return node;
            };
            System::ExecutionProgram program;program.root=parse(parse,object.parameters.at("root"));
            if (program.root.kind!=System::ExecutionKind::Sequence || program.root.children.empty())
                throw std::runtime_error("Authored HOW root requires a nonempty Sequence.");
            authoredExecution=std::move(program);
        });
    factories.add("providerBindings",{{{"bindings","array",true}},false},
        [&](const Model::ObjectDescriptor& object,Context&) {
            if (numericalBindingsDeclared) throw std::runtime_error("Multiple authored WHICH binding collections.");
            numericalBindingsDeclared=true;
            for (const auto& value:object.parameters.at("bindings")) {
                Model::Schema{{{"equation","string",true},{"method","string",true},{"inputs","array"},
                    {"occurrence","string"},{"parameters","object"}},false}.validate(value,"WHICH binding");
                System::NumericalBinding binding;
                binding.equation=value.at("equation").get<std::string>();binding.method=value.at("method").get<std::string>();
                binding.occurrence=value.value("occurrence",std::string());
                if (value.contains("inputs")) binding.inputs=value.at("inputs").get<std::vector<std::string>>();
                if (value.contains("parameters")) for (const auto& entry:value.at("parameters").items()) {
                    if (!entry.value().is_number()) throw std::runtime_error("WHICH parameters must be numeric: "+entry.key());
                    binding.parameters.emplace(entry.key(),entry.value().get<double>());
                }
                authoredNumerics.push_back(std::move(binding));
            }
        });
    factories.add("algorithmRegistry",{{{"Explicit","object"},{"SIMPLE","object"},
        {"PISO","object"},{"PIMPLE","object"}},false},
        [&](const Model::ObjectDescriptor& object,Context&) {
            composition.declared = true;
            if (!composition.algorithm.empty()) throw std::runtime_error("Multiple HOW algorithm selections.");
            int selected = 0;
            const auto select = [&](const char* name) {
                if (!object.parameters.contains(name)) return;
                if (++selected != 1) throw std::runtime_error(
                    "algorithms.yaml must select exactly one algorithm preset.");
                composition.algorithm = name;
                const auto& parameters = object.parameters.at(name);
                nativeAlgorithmParameters=parameters;
                if (!parameters.is_object()) throw std::runtime_error(
                    std::string("Algorithm '")+name+"' must be a mapping.");
                composition.outerCorrectors = parameters.value("outerCorrectors",1);
                composition.pressureCorrectors = parameters.value(
                    composition.algorithm == "PISO" ? "correctors" : "pressureCorrectors",1);
                composition.nonOrthogonalCorrectors = parameters.value(
                    "nonOrthogonalCorrectors",0);
                if (composition.outerCorrectors <= 0 || composition.pressureCorrectors <= 0
                    || composition.nonOrthogonalCorrectors < 0) throw std::runtime_error(
                    std::string("Algorithm '")+name+"' has invalid corrector counts.");
            };
            select("Explicit");
            select("SIMPLE");
            select("PISO");
            select("PIMPLE");
            if (selected == 0) throw std::runtime_error(
                "algorithms.yaml must select one registered algorithm preset.");
        });
    // Native object type → typed semantic slot. The sink name is the model's own
    // semantic role, never a disk layout path.
    static const std::map<std::string,std::string> typedRoles{
        {"Thermophysical","Thermophysical"},{"PhaseSystem","PhaseSystem"},
        {"PhaseChange","PhaseChange"},{"Turbulence","Turbulence"},
        {"ImmersedBoundary","ImmersedBoundary"},{"ILW","ILW"},
        {"Gravity","Gravity"},{"MRF","MRF"},{"WallHeat","WallHeat"},
        {"thermophysical","Thermophysical"},{"phaseChange","PhaseChange"},
        {"turbulence","Turbulence"},{"IBM","ImmersedBoundary"},
        {"gravity","Gravity"},{"wallHeat","WallHeat"}};
    const auto bindThermophysical=[&](const Model::ObjectDescriptor& o,Context& emit) {
        if (o.parameters.contains("thermoDynamics") || o.parameters.contains("properties")) {
            if(thermophysical)throw std::runtime_error("Duplicate thermophysical authority: "+thermophysical->source+" and "+o.name);
            thermophysical=parseThermophysical(o);
            composition.declared=true;
            composition.thermoDynamics=thermophysical->selection;
        } else if(o.type=="thermoDynamics")throw std::runtime_error("thermoDynamics object requires selection and properties.");
        emit("Thermophysical",o.parameters);
    };
    for(const auto& role:typedRoles) {
        if(role.second=="Thermophysical")factories.add(role.first,{{},true},bindThermophysical);
        else factories.add(role.first,{{},true},[role](const Model::ObjectDescriptor& o,Context& emit){emit(role.second,o.parameters);});
    }
    factories.add("thermoDynamics",{{},true},bindThermophysical);
    factories.add("multiPhase",{{},true},[](const Model::ObjectDescriptor& o,Context& emit) {
        auto parameters=o.parameters;
        if(!o.selection.empty()&&!parameters.contains("phaseSystem")) {
            parameters["phaseSystem"]=o.selection=="eulerian-eulerian"?"eulerianEulerian":o.selection;
        }
        emit("PhaseSystem",parameters);
    });
    // Equation expressions are deserialized generically. This solver only executes its registered presets.
    factories.add("BuiltinEquationSystem",{{{"preset","string",true}},false},[](const Model::ObjectDescriptor& o,Context&) {
        const auto preset=o.parameters.at("preset").get<std::string>();
        if(preset!="compressible"&&preset!="eulerianEulerian")throw std::runtime_error("Unknown equation preset: "+preset);
        if(!o.expression.empty())throw std::runtime_error("Custom equation expressions require an Equation compiler; built-in presets cannot ignore an expression.");
    });
    factories.add("equationSystem",{{{"preset","string",true}},false},[](const Model::ObjectDescriptor& o,Context&) {
        const auto preset=o.parameters.at("preset").get<std::string>();
        if(preset!="compressible"&&preset!="eulerianEulerian")throw std::runtime_error("Unknown equation preset: "+preset);
    });
    factories.add("File",{{{"file","string",true}},false},[&](const Model::ObjectDescriptor& o,Context&) {
        if(!std::filesystem::exists(std::filesystem::path(caseDir_)/o.parameters.at("file").get<std::string>()))throw std::runtime_error("Missing geometry: "+o.name);
    });
    factories.add("geometry",{{{"file","string",true}},false},[&](const Model::ObjectDescriptor& o,Context&) {
        if(!std::filesystem::exists(std::filesystem::path(caseDir_)/o.parameters.at("file").get<std::string>()))throw std::runtime_error("Missing geometry: "+o.name);
    });
    Context emit=assignRole;
    for(const auto& o:m.objects)factories.create(o,emit);
    for(const auto& role:roles) {
        const std::string& name=role.first;
        const P& value=role.second;
        if(name=="Thermophysical"){sections_.thermoDynamics=value;sections_.hasThermoDynamics=true;}
        else if(name=="PhaseSystem"){sections_.phaseSystem=value;sections_.hasPhaseSystem=true;}
        else if(name=="PhaseChange"){sections_.phaseChange=value;sections_.hasPhaseChange=true;}
        else if(name=="Turbulence"){sections_.turbulence=value;sections_.hasTurbulence=true;}
        else if(name=="ImmersedBoundary"){sections_.ibm=value;sections_.hasIbm=true;}
        else if(name=="ILW"){sections_.ilw=value;sections_.hasIlw=true;}
        else if(name=="Gravity"){sections_.gravity=value;sections_.hasGravity=true;}
        else if(name=="MRF"){sections_.mrf=value;sections_.hasMrf=true;}
        else if(name=="WallHeat"){sections_.wallHeat=value;sections_.hasWallHeat=true;}
        else throw std::runtime_error("Unknown built-in configuration role: "+name);
    }
    // field 的语义名字就是它的解码入口；不再合成 <startTime>/<name> 伪文档。
    for(const auto& f:m.fields) {
        // Derived views can have boundary laws without an independent initial value.
        if(f.initial.is_null()&&f.boundaries.empty())continue;
        Model::FieldDescriptor field=f;
        field.boundaries=expandBoundaryDefault(f.boundaries,caseDir_,meshFiles_);
        sections_.fields.push_back(std::move(field));
    }
    if (composition.declared) {
        if (!composition.stateDeclared || composition.solutionVariables.empty())
            throw std::runtime_error("Native composition requires explicit solution STATE (stateRegistry.use).");
        if (composition.algorithm.empty())
            throw std::runtime_error("Native composition requires explicit HOW (algorithmRegistry).");
        sections_.nativeAlgorithm=composition.algorithm;
        sections_.hasAlgorithm=true;
        // Existing linear/relaxation controls retain their decoder; family labels do not.
        P controls=sections_.algorithm.is_object()?sections_.algorithm:P::object();
        const auto oldName=controls.value("algorithm",std::string());
        if (!oldName.empty() && controls.contains(oldName) && controls.at(oldName).is_object())
            { const P oldControls=controls.at(oldName); controls.update(oldControls); }
        controls.erase("type");
        for (const auto* name:{"Explicit","SIMPLE","PISO","PIMPLE"}) controls.erase(name);
        controls.update(nativeAlgorithmParameters);
        if (controls.contains("correctors")) {
            controls["pressureCorrectors"]=controls.at("correctors");
            controls.erase("correctors");
        }
        P bound=P::object();
        bound["algorithm"]=composition.algorithm;
        if (controls.contains("linearSolvers")) {
            bound["linearSolvers"]=controls.at("linearSolvers");
            controls.erase("linearSolvers");
        }
        controls.erase("algorithm");
        bound[composition.algorithm]=std::move(controls);
        sections_.algorithm=std::move(bound);
    }
    if (authoredExecution && (!numericalBindingsDeclared || authoredNumerics.empty()))
        throw std::runtime_error("Authored HOW requires an explicit nonempty providerBindings WHICH collection; preset bindings are not silently inherited.");
    auto result=decodeNativeCase();if(!result)throw std::runtime_error("Cannot bind model "+m.name);
    result->caseName=m.name;
    if(result->createMesh)result->meshParameterFile=meshGenerator_;
    result->modelDescription=std::make_shared<const Model::Description>(m);
    if (!composition.declared) {
        // The old adapter emits an explicit tuple once. Runtime never reads its label.
        composition.stateDeclared=true;
        composition.solutionVariables={"rho","rhoU","rhoE"};
        composition.stateSelectionOrigin="compatibility "+legacyFlowLabel_;
        composition.equations={"Continuity","Momentum","Energy"};
        composition.algorithm=result->pressureCouplingDeclared
            ? FDM::toString(result->solver.pressure.coupling.preset) : "Explicit";
        composition.compatibilityPressureConstraint=legacyFlowLabel_=="pressureBase";
    }
    composition.outerCorrectors=result->solver.pressure.coupling.outerCorrectors;
    composition.pressureCorrectors=result->solver.pressure.coupling.pressureCorrectors;
    composition.nonOrthogonalCorrectors=result->solver.pressure.coupling.nonOrthogonalCorrectors;
    const bool fluid=std::find(composition.equations.begin(),composition.equations.end(),"Momentum")!=composition.equations.end();
    if (!thermophysical && !composition.declared && fluid && (!result->multiPhaseEnabled || Physics::Multiphase::isLevelSetType(result->multiPhase.type))) {
        ThermophysicalContract t;t.selection={"perfectGas","hConst","const",0.};
        t.source="legacy single-fluid normalization";
        t.gamma=result->solver.numerics.idealGasGamma;t.gasConstant=result->solver.numerics.idealGasConstant;
        t.dynamicViscosity=result->solver.numerics.dynamicViscosity;t.prandtl=result->solver.numerics.prandtl;
        t.provenance={{"gamma","historical NumericsConfig default 1.4"},{"R","historical NumericsConfig default 287.05"},
            {"mu",m.numerics.contains("transport")&&m.numerics.at("transport").contains("mu")?"legacy numerics.transport.mu":"historical mu default 0"},
            {"Pr",m.numerics.contains("transport")&&m.numerics.at("transport").contains("Pr")?"legacy numerics.transport.Pr":"historical Pr default 0.72"}};
        t.validate();thermophysical=t;
    }
    if(thermophysical) {
        if(thermophysical->native && m.numerics.contains("transport")) {
            const auto& legacy=m.numerics.at("transport");
            for(const auto* name:{"mu","Pr"})if(legacy.contains(name)) {
                const double v=legacy.at(name).get<double>();
                const double expected=std::string(name)=="mu"?thermophysical->dynamicViscosity:thermophysical->prandtl;
                if(!std::isfinite(v)||v!=expected)throw std::runtime_error("Native/legacy thermophysical parameter conflict: "+std::string(name));
                thermophysical->provenance[name]+="; equal legacy numerics.transport."+std::string(name);
            }
        }
        auto& n=result->solver.numerics;
        if(thermophysical->selection.equationOfState=="perfectGas") {
            n.idealGasGamma=thermophysical->gamma;n.idealGasConstant=thermophysical->gasConstant;
        }
        n.dynamicViscosity=thermophysical->dynamicViscosity;n.prandtl=thermophysical->prandtl;
        result->solver.boundaries.thermalDynamicViscosity=n.dynamicViscosity;
        result->solver.boundaries.thermalPrandtl=n.prandtl;
        result->solver.turbulence.laminarDynamicViscosity=n.dynamicViscosity;
        const bool selectedDiffusion=n.recipes.diffusion.has_value()
            || (m.numerics.contains("diffusion") && m.numerics.at("diffusion").contains("default"));
        // Native operator selection drives the old enable projection; transport
        // coefficients no longer have to masquerade as numerical configuration.
        const bool legacyDiffusionDisabled=m.numerics.contains("transport")
            && m.numerics.at("transport").contains("enable")
            && !m.numerics.at("transport").at("enable").get<bool>();
        if (thermophysical->native && selectedDiffusion && !legacyDiffusionDisabled) n.viscousEnabled=true;
        if (legacyDiffusionDisabled)
            thermophysical->provenance["diffusion activation"]="legacy numerics.transport.enable=false; declared scheme remains inactive";
        if (thermophysical->native && n.viscousEnabled) {
            if (thermophysical->selection.transport.empty())
                throw std::runtime_error(thermophysical->source+": diffusion requires a selected transport provider.");
            if (thermophysical->selection.equationOfState=="perfectGas" && thermophysical->dynamicViscosity==0.)
                throw std::runtime_error(thermophysical->source+": current viscous heat-flux provider requires positive conductivity; mu=0 with diffusion is unsupported.");
        }
        result->solver.thermophysical=std::make_shared<const ThermophysicalContract>(*thermophysical);
        result->solver.boundaries.thermophysical=result->solver.thermophysical;
        result->solver.validateThermophysicalProjection();
    }
    result->composition=std::move(composition);
    result->authoredExecution=std::move(authoredExecution);
    result->authoredNumerics=std::move(authoredNumerics);
    return *result;
}
}
