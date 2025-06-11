if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG fdaac01bcc52888994f7afd029dcc045dd408484 # v2.5.0
)
FetchContent_MakeAvailable(libigl)
