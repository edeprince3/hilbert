foreach(dir 
    liblbfgs
    libsdp
)
    add_subdirectory(${dir})
endforeach()
