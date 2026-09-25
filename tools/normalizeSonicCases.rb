#!/usr/bin/env ruby
# frozen_string_literal: true

# Convert the repository's transitional native cases to the canonical SonicFile layout.
# This script deliberately emits the user-facing SonicFile syntax, including
# `patch: fixedValue: value`; SF_serialization normalizes that convenience syntax before YAML parsing.
require 'yaml'
require 'pathname'

ROOT = Pathname.new(__dir__).join('..').realpath
TYPE_MAP = {
  'BuiltinEquationSystem' => 'equationSystem', 'ImmersedBoundary' => 'IBM',
  'PhaseSystem' => 'multiPhase', 'Thermophysical' => 'thermophysical',
  'PhaseChange' => 'phaseChange', 'Turbulence' => 'turbulence',
  'Gravity' => 'gravity', 'WallHeat' => 'wallHeat', 'File' => 'geometry'
}.freeze
PREFERRED_NAME = {
  'PhaseSystem' => 'multiPhase', 'Thermophysical' => 'thermophysical',
  'PhaseChange' => 'phaseChange', 'Turbulence' => 'turbulence',
  'Gravity' => 'gravity', 'MRF' => 'MRF', 'WallHeat' => 'wallHeat'
}.freeze

def load_yaml(path)
  YAML.load_file(path) || {}
end

def body_of(data)
  return data unless data.is_a?(Hash)
  root = data.dup
  if root['SonicFile'].is_a?(Hash)
    root.delete('SonicFile')
  else
    root.delete('object')
    root.delete('type')
  end
  root
end

def plain_string?(value)
  value.match?(/\A[A-Za-z_][A-Za-z0-9_\.\/+-]*\z/) && !%w[true false null yes no on off].include?(value.downcase)
end

def scalar(value)
  case value
  when NilClass then 'null'
  when TrueClass then 'true'
  when FalseClass then 'false'
  when Integer then value.to_s
  when Float
    raise "nonfinite value #{value}" unless value.finite?
    value.to_s
  when String
    # Psych 3 writes a scalar such as "1e-08" as a complete YAML document ("--- 1e-08").
    # Keep numeric text bare so YAML restores its numeric type, and quote only true text strings.
    return value if value.match?(/\A[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?\z/)
    plain_string?(value) ? value : value.inspect
  else raise "unsupported scalar #{value.class}"
  end
end

def inline(value)
  case value
  when Array then '[' + value.map { |item| inline(item) }.join(', ') + ']'
  when Hash
    '{' + value.map { |key, item| "#{key}: #{inline(item)}" }.join(', ') + '}'
  else scalar(value)
  end
end

def emit_mapping(io, body, indent = 0, compact_conditions = false)
  raise "SonicFile body must be a mapping" unless body.is_a?(Hash)
  pad = ' ' * indent
  body.each do |key, value|
    if value.is_a?(Hash)
      if value.empty?
        io.puts "#{pad}#{key}: {}"
      elsif compact_conditions && value.length == 1 && !value.values.first.is_a?(Hash)
        condition, prescribed = value.first
        suffix = prescribed.nil? ? '' : ": #{inline(prescribed)}"
        io.puts "#{pad}#{key}: #{condition}#{suffix}"
      else
        io.puts "#{pad}#{key}:"
        emit_mapping(io, value, indent + 2, compact_conditions)
      end
    elsif value.is_a?(Array)
      io.puts "#{pad}#{key}: #{inline(value)}"
    else
      io.puts "#{pad}#{key}: #{scalar(value)}"
    end
  end
end

def write_sonic(path, object, type, body, comment = nil, compact_conditions: false)
  path.dirname.mkpath
  path.open('w') do |io|
    io.puts "# #{comment}" if comment
    io.puts 'SonicFile:'
    io.puts "  object: #{object}"
    io.puts "  type: #{type}"
    emit_mapping(io, body, 0, compact_conditions)
  end
end

def unique_name(seed, used)
  base = seed.gsub(/[^A-Za-z0-9_-]/, '_')
  base = 'module' if base.empty?
  name = base
  index = 2
  while used.include?(name)
    name = "#{base}#{index}"
    index += 1
  end
  used << name
  name
end

Dir.glob(ROOT.join('test/**/case.yaml').to_s).sort.each do |case_file|
  directory = Pathname.new(case_file).dirname
  next if directory.to_s.include?('/ioRegistry/')
  legacy_case = load_yaml(case_file)
  write_sonic(directory.join('case.yaml'), 'case', 'registry', body_of(legacy_case), 'Case identity.')

  mesh = body_of(load_yaml(directory.join('mesh/mesh.yaml')))
  write_sonic(directory.join('mesh/mesh.yaml'), 'mesh', 'registry', mesh, 'Mesh registry.')

  fields = body_of(load_yaml(directory.join('fields/fields.yaml')))
  initial = body_of(load_yaml(directory.join('fields/internalField.yaml')))
  boundaries = body_of(load_yaml(directory.join('fields/boundaries.yaml')))
  write_sonic(directory.join('fields/fields.yaml'), 'fields', 'registry', fields, 'Field registry.')
  write_sonic(directory.join('fields/internalField.yaml'), 'fields', 'internalField', initial, 'Initial field values.')
  write_sonic(directory.join('fields/boundaries.yaml'), 'fields', 'boundary', boundaries,
              'Boundary conditions. SonicFile accepts patch: condition: value.', compact_conditions: true)

  solvers = body_of(load_yaml(directory.join('solvers/solvers.yaml')))
  models = body_of(load_yaml(directory.join('models/models.yaml')))
  solver_registry = {
    'runtime' => { 'type' => 'runtime', 'file' => 'solvers/runtime.yaml' },
    'numerics' => { 'type' => 'numerics', 'file' => 'solvers/numerics.yaml' },
    'algorithm' => { 'type' => 'algorithm', 'file' => 'solvers/algorithm.yaml' }
  }
  write_sonic(directory.join('solvers/runtime.yaml'), 'solver', 'runtime', solvers.fetch('runtime', {}), 'Runtime configuration.')
  write_sonic(directory.join('solvers/numerics.yaml'), 'solver', 'numerics', solvers.fetch('numerics', {}), 'Numerical methods.')
  write_sonic(directory.join('solvers/algorithm.yaml'), 'solver', 'algorithm', solvers.fetch('solver', {}), 'Solver algorithm.')

  equation = solvers['equation'] || models.dig('equations', 'flow')
  if equation
    parameters = equation['parameters'] || {}
    write_sonic(directory.join('models/equations.yaml'), 'solver', 'equationSystem', parameters, 'Flow equation system.')
    solver_registry['flow'] = { 'type' => 'equationSystem', 'file' => 'models/equations.yaml' }
  end

  %w[IBM ILW].each do |name|
    source = directory.join("solvers/#{name}.yaml")
    next unless source.file?
    module_body = body_of(load_yaml(source))
    write_sonic(directory.join("models/#{name}.yaml"), 'solver', name, module_body, "Built-in #{name} module.")
    # Built-ins use the IO registry's declared default path. Custom entries always carry file.
    solver_registry[name] = { 'type' => name }
    source.delete
  end
  write_sonic(directory.join('solvers/solvers.yaml'), 'solver', 'registry', solver_registry, 'Solver registry.')

  model_registry = {}
  used_names = []
  %w[closures constraints geometry].each do |category|
    entries = models[category]
    next unless entries.is_a?(Hash)
    entries.each do |legacy_name, descriptor|
      next unless descriptor.is_a?(Hash) && descriptor['type']
      type = descriptor['type']
      semantic = TYPE_MAP.fetch(type, type)
      seed = PREFERRED_NAME.fetch(type, legacy_name)
      name = unique_name(seed, used_names)
      parameters = descriptor['parameters'] || {}
      body = descriptor['expression'] ? { 'expression' => descriptor['expression'], 'parameters' => parameters } : parameters
      file = "models/#{name}.yaml"
      write_sonic(directory.join(file), 'models', semantic, body, "#{name} model module.")
      model_registry[name] = { 'type' => semantic, 'file' => file }
    end
  end
  write_sonic(directory.join('models/models.yaml'), 'models', 'registry', model_registry, 'Physics model registry.')
end
