# bitflip change
#export PETSC_DIR=/import/freenas-m-03-geodynamics/jhayek/petsc-3.12.5

#export PETSC_DIR=/import/freenas-m-02-seismology/jhayek/TEAR/Project_01_Petsc/petsc
#export PETSC_ARCH=arch-linux-c-debug
#export LD_LIBRARY_PATH=${PETSC_DIR}/${PETSC_ARCH}/lib

cd ../se2wave
make clean
make all

now=$(date)
echo ">> $now"

export STORAGE_VAR_NAME=/home/nico/Documents/TEAR/Codes_TEAR/se2dr/Runs/Tilt20_50x50_25delta_Force
export SE2DR_RUN="/home/nico/Documents/TEAR/Codes_TEAR/se2dr/se2wave/se2dr.app"
mkdir -p $STORAGE_VAR_NAME

#(cd $STORAGE_VAR_NAME && $SE2DR_RUN -mx 160 -my 160 -tmax 4 -bdegree 1 -delta_cell_factor 1.001 -of 10)

#(cd $STORAGE_VAR_NAME && $SE2DR_RUN -mx 200 -my 200 -tmax 4 -bdegree 1 -delta_cell_factor 1.001 -of 10)
#(cd $STORAGE_VAR_NAME && $SE2DR_RUN -mx 400 -my 400 -tmax 4 -bdegree 1 -delta_cell_factor 1.001 -of 10)
(cd $STORAGE_VAR_NAME && $SE2DR_RUN -mx 400 -my 400 -tmax 4 -bdegree 3 -delta 25 -of 10)

# Debug~~
#valgrind --tool=memcheck --leak-check=full ./se2dr.app -mx 200 -my 200 -tmax 4.0 -bdegree 2 -delta_cell_factor 1.001 -of 10
