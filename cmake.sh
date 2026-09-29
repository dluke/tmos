cmake -S code -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER=/usr/bin/c++ \
    -DPython3_ROOT_DIR="$HOME/.pyenv/versions/3.11.16" \
    -DPython3_EXECUTABLE="$HOME/.pyenv/versions/3.11.16/bin/python3.11" \
    -DPython3_LIBRARY="$HOME/.pyenv/versions/3.11.16/lib/libpython3.11.so" \
    -DUSE_VTK=ON \
    -DVTK_DIR="$HOME/lib/vtk-build/lib/cmake/vtk-9.3" \
    -DBOOST_ROOT=$HOME/lib/boost_1_67_0 \
    -DBoost_LIBRARY_DIR="$HOME/lib/boost_1_67_0/stage/lib" \
