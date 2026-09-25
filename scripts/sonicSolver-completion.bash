# sonicSolver bash completion.
# Install/source this file after building with ./sonicMake.sh.

_sonicSolver_complete() {
    local cur prev command
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    command="${COMP_WORDS[1]}"

    local commands="run check explain doctor models why explain-model recipes init cleanCase postProcess -cleanResult -postProcessing"
    local options="--help -h --initial-output --steps"
    local models="ghostCell peskinOriginal dfmExplicitSelfPropelled dfmFractionalStepSelfPropelled dfmFractionalStepPrescribed dfmImplicitPrescribed dfmImplicitSelfPropelled dfmAugmentedLagrangian velocityForcingFTS velocityForcingBP"
    local recipes="compressibleSingleFluid compressibleFixedImmersedBody compressibleMovingBody selfPropelledBody levelSetSurfaceTension"

    if [[ "$prev" == "--with" ]]; then
        COMPREPLY=( $(compgen -W "mpi" -- "$cur") )
        return 0
    fi
    if [[ "$prev" == "--recipe" ]]; then
        COMPREPLY=( $(compgen -W "$recipes" -- "$cur") )
        return 0
    fi

    case "$command" in
        models)
            COMPREPLY=( $(compgen -W "ibm --with" -- "$cur") )
            ;;
        why|explain-model)
            COMPREPLY=( $(compgen -W "$models" -- "$cur") )
            ;;
        recipes)
            COMPREPLY=()
            ;;
        init)
            if [[ "$cur" == -* ]]; then
                COMPREPLY=( $(compgen -W "--recipe" -- "$cur") )
            elif (( COMP_CWORD == 2 )); then
                COMPREPLY=( $(compgen -d -- "$cur") )
            else
                COMPREPLY=( $(compgen -W "--recipe" -- "$cur") )
            fi
            ;;
        run|check|explain|doctor|cleanCase|postProcess|-cleanResult|-postProcessing|--initial-output|--steps)
            if [[ "$cur" == -* ]]; then
                COMPREPLY=( $(compgen -W "$options" -- "$cur") )
            else
                COMPREPLY=( $(compgen -d -- "$cur") )
            fi
            ;;
        *)
            if [[ "$cur" == -* ]]; then
                COMPREPLY=( $(compgen -W "$options" -- "$cur") )
            else
                COMPREPLY=( $(compgen -W "$commands" -- "$cur") $(compgen -d -- "$cur") )
            fi
            ;;
    esac
}

complete -F _sonicSolver_complete sonicSolver
