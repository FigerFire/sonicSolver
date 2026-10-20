# sonicSolver bash completion.
# Install/source this file after building with ./sonicMake.sh.

_sonicSolver_complete() {
    local cur prev command
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    command="${COMP_WORDS[1]}"

    local commands="run check explain doctor models why explainModel recipes init cleanResult postProcessing initialOutput steps help"
    local options="--initial-output --steps"
    local models="ghostCell peskinOriginal dfmExplicitSelfPropelled dfmFractionalStepSelfPropelled dfmFractionalStepPrescribed dfmImplicitPrescribed dfmImplicitSelfPropelled dfmAugmentedLagrangian velocityForcingFTS velocityForcingBP"
    local recipes="compressibleSingleFluid compressibleFixedImmersedBody compressibleMovingBody selfPropelledBody levelSetSurfaceTension"

    COMPREPLY=()
    if (( COMP_CWORD == 1 )); then
        COMPREPLY=( $(compgen -W "$commands" -- "$cur") )
        return 0
    fi

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
        why|explainModel)
            COMPREPLY=( $(compgen -W "$models" -- "$cur") )
            ;;
        steps)
            if (( COMP_CWORD > 2 )); then
                COMPREPLY=( $(compgen -d -- "$cur") )
            fi
            ;;
        help|recipes)
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
        run)
            if [[ "$prev" == "--steps" ]]; then
                COMPREPLY=()
            elif [[ "$cur" == -* ]]; then
                COMPREPLY=( $(compgen -W "$options" -- "$cur") )
            else
                COMPREPLY=( $(compgen -d -- "$cur") )
            fi
            ;;
        check|explain|doctor|cleanResult|postProcessing|initialOutput)
            if [[ "$cur" == -* ]]; then
                COMPREPLY=()
            else
                COMPREPLY=( $(compgen -d -- "$cur") )
            fi
            ;;
        *)
            if [[ "$cur" == -* ]]; then
                COMPREPLY=()
            else
                COMPREPLY=( $(compgen -d -- "$cur") )
            fi
            ;;
    esac
}

complete -o filenames -F _sonicSolver_complete sonicSolver
