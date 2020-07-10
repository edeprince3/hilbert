# hilbert

**Installation**

1. install psi4. i like to install the "nightly build" from source (see http://vergil.chemistry.gatech.edu/nu-psicode/install.html )

2. download hilbert:

    git clone git@github.com:edeprince3/hilbert.git
  
3. configure and compile hilbert:

    psi4 --plugin-compile > do-configure
    
    (add any additional cmake flags to do-configure you'd like)
    
    chmod +x do-configure
    
    ./do-configure
    
    make
  
4. test:

    cd tests
    make
