"""Targeted STATE/math/scalar/provider checks; every working directory is below test/."""
from pathlib import Path
import json,os,subprocess
base=Path(__file__).resolve().parent
root=base.parents[1]
os.chdir(base)
env={**os.environ,'DEVELOPER_DIR':'/Library/Developer/CommandLineTools','TMPDIR':str(base/'evidence')}
selected={
'mathSignatures','typedWhat','mathNumericalBinding','formulaExecutionContract',
'compiledResourceComposition','coreSinglePathClosure','couplingContributionsLifecycle','derivedStateResourceViews','standaloneScalarExecution','coupledScalarExecution','multiScalarExecution','nonlinearScalarExecution',
'runtimeBindingLifetime','explicitStageMathematics','mixedTemporalComposition','scopeAndStateContracts',
'formulationArchitecture','methodObjectCompilation','providerCatalogResolution','termRecipeAuthority',
'fluidScalarExecution','sharedFluxTransport','stageSSTTransport','stageCoefficientContracts','stageCoefficientEvidence',
'pisoArchitecture','fixedTimePressureContract','nativeEulerianExecution','nativeEulerianTurbulenceContract',
'nativeThermalExecution','nativeLevelSetExecution','nativeMultiphaseCapability','nativeTurbulenceExecution',
'rasOnceStepOwnership','lesPhysicalDensity','thermophysicalAuthority','builtinCompositionClosure',
'programDoesNotOwnRuntimeState','pressureTermProviders','numericalFluxContract','nativeImmersedContributions','viscousImmersedWall','surfaceSchurProjection'
}
inventory=json.loads(subprocess.check_output(['ctest','--test-dir',str(root/'build'),'--show-only=json-v1'],env=env,text=True))
q=lambda v:json.dumps(str(v));lines=[];found=[]
for item in inventory['tests']:
 if item['name'] not in selected:continue
 name=item['name'];found.append(name)
 command=list(item['command'])
 if name=='nonlinearScalarExecution':command[-1]=str(base/'cases/nonlinear')
 cwd=base/'evidence/unit'/name;cwd.mkdir(parents=True,exist_ok=True)
 lines+=['add_test('+q(name)+' '+' '.join(q(c) for c in command)+')','set_tests_properties('+q(name)+' PROPERTIES WORKING_DIRECTORY '+q(cwd)+' TIMEOUT '+('180' if name in {'builtinCompositionClosure','pressureTermProviders'} else '15')+')']
assert set(found)==selected,selected-set(found)
folder=base/'ctest';folder.mkdir(exist_ok=True);(folder/'CTestTestfile.cmake').write_text('\n'.join(lines))
(base/'evidence/test-selection.json').write_text(json.dumps(found,indent=2))
with (base/'evidence/ctest-final.log').open('w') as log:
 result=subprocess.run(['ctest','--test-dir',str(folder),'--output-on-failure','--parallel','1'],env=env,stdout=log,stderr=subprocess.STDOUT)
print((base/'evidence/ctest-final.log').read_text())
raise SystemExit(result.returncode)
