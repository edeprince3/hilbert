# hilbert

**Installation**

1. install psi4. i like to install the "nightly build" from source using an appropriate conda environment (see http://vergil.chemistry.gatech.edu/nu-psicode/install.html ).

2. download hilbert:

        git clone git@github.com:edeprince3/hilbert.git
  
3. use psi4 to generate an appropriate cmake configure line:

        conda activate p4dev
        psi4 --plugin-compile 
    
    add any additional cmake flags to the psi4-generated cmake line that you'd like. the orbital optimizer generates quite a few temporary files, so consider compiling in a separate directory. in that case, you should also set the install directory, i.e., 
    
        cmake ... -Bobjdir -DCMAKE_INSTALL_PREFIX=path_to_hilbert/hilbert  
        cd objdir    
        make -j N   
        make install
        cd ..
  
4. test:

        cd tests
        make
