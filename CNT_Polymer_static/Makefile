objects = data_structure.o main.o readin.o Shape_function.o Sparse_matrix.o Tcell.o Tcomposite.o\
               TFE_grid.o TMKL_solver.o
CNTcomposite_simulation : $(objects)
	icc -o CNTcomposite_simulation $(objects) -mkl=parallel -qopenmp
data_structure.o : data_structure.cpp
	icc  -c -O3 -xHost data_structure.cpp
main.o : main.cpp Tcomposite.h
	icc  -c -O3 -xHost  main.cpp
readin.o : readin.cpp data_structure.h
	icc  -c -O3 -xHost  readin.cpp
Shape_function.o : Shape_function.cpp data_structure.h
	icc  -c -O3 -xHost  Shape_function.cpp
Sparse_matrix.o : Sparse_matrix.cpp
	icc  -c -O3 -xHost Sparse_matrix.cpp
Tcell.o : Tcell.cpp data_structure.h Shape_function.h Gauss_weight.h
	icc  -c -O3 -xHost Tcell.cpp
Tcomposite.o : Tcomposite.h
	icc  -c -O3 -xHost Tcomposite.cpp -qopenmp -lpthread
TFE_grid.o : TFE_grid.cpp data_structure.h Tcell.h readin.h TMKL_solver.h Sparse_matrix.h Shape_function.h
	icc  -c -O3 -xHost  TFE_grid.cpp
TMKL_solver.o : TMKL_solver.cpp 
.PHONY : clean
clean : 
	rm CNTcomposite_simulation $(objects)
