
#
# cleanup
#

tar_file:	all_clean
		tar cvf op2.tar \
		    jac1/*.h jac1/*.cpp jac1/*.cu jac1/Makefile \
		    jac2/*.h jac2/*.cpp jac2/*.cu jac2/Makefile \
		    airfoil/*.h airfoil/*.cpp airfoil/*.cu airfoil/Makefile \
		    common/*.h common/*.cpp common/*.cu common/*.m \
                    doc/*.tex doc/*.pdf \
                    README Makefile

all_clean:	jac1_clean jac2_clean airfoil_clean


jac1_clean:	
		cd jac1; make clean;

jac2_clean:	
		cd jac2; make clean;

airfoil_clean:	
		cd airfoil; make clean;


