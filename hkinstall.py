from hkpilot.utils import fancylogger
from hkpilot.utils.cmake import CMake

logger = fancylogger.getLogger(__name__)

class FLOWER(CMake):

    def __init__(self, path):
        super().__init__(path)

        self._package_name = "FLOWER"
        self._cmakelist_path = "."
        self._cmake_options = {
            "CMAKE_CXX_STANDARD": "14",
            "HEMI_CUDA_DISABLE": "ONE"
        }
