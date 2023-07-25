#!/bin/bash

# Input parsing --------------------------------------------------------------

function usage {
printf "Usage: ./run_study.sh [OPTION]...
Run parameter sweep of fun3d runs for agard wing.

\t-g, --get_grid   download agard 445.6 grid from fun3d website
\t-p, --prepare    prepare new files/folders
\t-r, --run        run complete parameter sweep (1=startup, 2=limit-cycle)
"
}

while getopts 'gpr:' OPTION; do
    case "$OPTION" in
        g)
            get_grid=true
            ;;
        p)
            prepare=true
            ;;
        r)
            run=true
            run_type=$OPTARG
            ;;
        ?)
            usage
            exit 1
            ;;
    esac
done
shift "$(($OPTIND - 1))"

# No input options - report usage
if (( $OPTIND == 1 )); then
    usage
    exit 0
fi

# User inputs ----------------------------------------------------------------
# run_indexes=(0 1 2 3 4 5)  # indexes of mach numbers to submit for simulation
run_indexes=(3)  # indexes of mach numbers to submit for simulation

# experimental
machs=(0.50 0.68 0.90 0.96 1.07 1.14)  # -, mach number
# speed_index_bases=(0.450 0.430 0.345 0.190 0.480 0.685)  # -, speed index
speed_index_bases=(0.450 0.430 0.345 0.240 0.480 0.685)  # -, speed index
densities=(0.000830 0.000404 0.000193 0.000123 0.000107 0.000152)  # slug/ft^3, density
mass_ratios=(33.465 68.753 143.920 225.820 259.590 182.740)  # -, mass ratio (structural/fluid)
dynamic_pressures=(135.3 123.5 79.5 24.1 153.9 313.5)  # psf, dynamic pressure
velocities=(570.9 782.0 907.7 626.2 1696.1 2030.9)  # ft/s, freestream velocity
L_star=0.9165  # ft, root semi-chord
L=L_star  # -, grid is in ft
omega=239.3  # rad/s, first uncoupled torsional natural frequency

# physics
time_step_nondim=0.3
reynolds=0.0  # -, fun3d default for inviscid
temperature=491.4  # R, fun3d default for inviscid (converted to Rankine)

# simulation
cores_per_run=8
speed_index_increments=(0.9 1.0 1.1)
run_name="wing-445.6"
step1="steady"  # run the steady-state flow case in fun3d
step2="modes"  # map the mode shapes from the FEA mesh to the CFD grid
step3="unsteady"  # run the unsteady, coupled fsi flutter case in fun3d

main_path=$(pwd)

# GET AGARD GRID -------------------------------------------------------------
if [ "$get_grid" = true ]; then
    wget https://fun3d.larc.nasa.gov/FUN3D_v13.4_session16.tar.gz
    tar -xvf FUN3D_v13.4_session16.tar.gz
    tar -xvf session16/flow_modal_aeroelasticity.tar.gz
    cp -r flow_modal_aeroelasticity/Grid grid
    rm -rf flow_modal_aeroelasticity session16 FUN3D_v13.4_session16.tar.gz
    ls
    exit 0
fi

# STEP 2 - MAP MODE SHAPES ---------------------------------------------------
modes_path=$main_path/results/$step2
grid_path=$main_path/results/grid
if [ "$prepare" = true ]; then

    # Create main results folder
    mkdir -pv $modes_path

    # Create grid folder
    mkdir -pv $grid_path

    # Copy grid
    cp -r ./grid/* $grid_path

    # Copy mode-mapping files/folder
    cp -r ./template/$step2/* $modes_path

fi
if [ "$run" = true ] && [ "$run_type" == "2" ]; then  # map mode shapes
    hm=$(pwd)
    cd $modes_path
    pwd

    # Link to massoud template file from steady solution
    run_ind=${run_indexes[0]}
    speed_ind=$(echo "scale=8;${speed_index_bases[$run_ind]}*${speed_index_increments[$run_ind]}" | bc -l)  # adjust speed index by small increment
    ln -s $main_path/results/mach-${machs[$run_ind]}/speed_index-${speed_ind}/steady/*massoud_body1.dat .

    # Compile mode shape mapping code
    gfortran -fdefault-real-8 Mode.f

    # Execute mode shape mapping code
    ./a.out > fun3d.out

    cd $hm
fi

for run_index in ${run_indexes[@]}; do  # loop through mach numbers
    mach=${machs[run_index]}
    echo $mach
    speed_index_base=${speed_index_bases[run_index]}

    density=${densities[run_index]}
    mass_ratio=${mass_ratios[run_index]}

    for speed_index_increment in ${speed_index_increments[@]}; do  # loop through speed indexes
        speed_index=$(echo "scale=8;$speed_index_base*$speed_index_increment" | bc -l)  # adjust speed index by small increment
        echo "  $speed_index"

        steady_path=$main_path/results/mach-${mach}/speed_index-${speed_index}/steady
        unsteady_path=$main_path/results/mach-${mach}/speed_index-${speed_index}/unsteady
        if [ "$prepare" = true ]; then

            # Calculate relevant parameters
            velocity=$(echo "scale=8;$speed_index*$L_star*$omega*sqrt($mass_ratio)" | bc -l)
            dynamic_pressure=$(echo "scale=8;0.5*$density*$velocity^2" | bc -l)
            # reynolds=$(echo "scale=8;$density*$velocity/$mu*($L_star/$L)" | bc -l)
            # speed_of_sound=$(echo "scale=8;$velocity/$mach" | bc -l)
            # pressure=$(echo "scale=8;$speed_of_sound^2*$density/$gamma" | bc -l)
            # temperature=$(echo "scale=8;$pressure/($density*$R)" | bc -l)

            # STEP 1 - STEADY ------------------------------------------------

            # Create steady run folder
            mkdir -pv $steady_path

            # Copy steady files/folder
            cp -r ./template/$step1/* $steady_path

            # Set mach number
            sed -i "s/mach_number_value/$mach/" $steady_path/fun3d.nml

            # Set reynolds number (dimensionalizes time)
            sed -i "s/reynolds_number_value/$reynolds/" $steady_path/fun3d.nml

            # Set temperature
            sed -i "s/temperature_value/$temperature/" $steady_path/fun3d.nml

            # STEP 3 - UNSTEADY ----------------------------------------------

            # Create unsteady run folder
            mkdir -pv $unsteady_path
            
            # Copy unsteady files/folder
            cp -r ./template/$step3/* $unsteady_path
            
            # Set mach number
            sed -i "s/mach_number_value/$mach/" $unsteady_path/fun3d.nml

            # Set reynolds number (dimensionalizes time)
            sed -i "s/reynolds_number_value/$reynolds/" $unsteady_path/fun3d.nml

            # Set temperature
            sed -i "s/temperature_value/$temperature/" $unsteady_path/fun3d.nml

            # Set time step (nondimensional)
            sed -i "s/time_step_nondim_value/$time_step_nondim/" $unsteady_path/fun3d.nml

            # Set freestream velocity (dimensionalizes time)
            sed -i "s/uinf_value/$velocity/" $unsteady_path/moving_body.input

            # Set dynamic pressure (dimensionalizes/scales modal forces)
            sed -i "s/qinf_value/$dynamic_pressure/" $unsteady_path/moving_body.input
            
        fi
        if [ "$run" = true ]; then

            hm=$(pwd)
            module load fun3d/13.7-tecio

            # STEP 1 - STEADY ------------------------------------------------

            if [ "$run_type" == "1" ]; then  # run steady solution
                cd $steady_path

                # Link to grid files
                ln -s $grid_path/*$run_name* ./

                # Submit fun3d run
                /usr/bin/time -a --output=fun3d.out --format='%e' mpiexec -np $cores_per_run nodet_mpi --write_massoud_file > fun3d.out &

            # STEP 3 - UNSTEADY ----------------------------------------------

            elif [ "$run_type" == "3" ]; then  # run unsteady solution
                cd $unsteady_path

                # Link to grid files
                ln -s $grid_path/*${run_name}* ./

                # Link to restart .flow file
                ln -s $steady_path/${run_name}.flow ./

                # Link to mapped mode shape files
                ln -s $modes_path/${run_name}_body1_mode* .

                # Submit fun3d run
                /usr/bin/time -a --output=fun3d.out --format='%e' mpiexec -np $cores_per_run nodet_mpi --aeroelastic_internal > fun3d.out &

            fi

            cd $hm

        fi
    done
done
