# hilbert

**Installation**

1. install psi4. i like to install the "nightly build" from source (see http://vergil.chemistry.gatech.edu/nu-psicode/install.html )

2. download hilbert:

    git clone git@github.com:edeprince3/hilbert.git
  
3. configure and compile hilbert:

    psi4 --plugin-compile > do-configure
    
    add any additional cmake flags to do-configure you'd like. the orbital optimizer generates some temporary files, so consider compiling in a separate directory. in that case, you should also set the install directory, i.e., 
    
        cmake ... -Bobjdir -DCMAKE_INSTALL_PREFIX=path_to_hilbert/hilbert
    
    chmod +x do-configure
    
    ./do-configure
    
    cd objdir
    
    make
    
    make install
  
4. test:

    cd tests
    make
