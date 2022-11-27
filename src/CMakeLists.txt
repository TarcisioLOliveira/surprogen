
set(SOURCES
    main.cpp
)

add_executable(${PROJECT_NAME} ${SOURCES})

add_subdirectory(solver)
add_subdirectory(engine2D)

target_link_libraries(${PROJECT_NAME}
    engine2D
    solver 
    Python3::Python
    Python3::Module
    Python3::NumPy
    ${YAML_CPP_LIBRARIES}
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)
endif()

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION .
                                LIBRARY DESTINATION .
                                ARCHIVE DESTINATION .)
