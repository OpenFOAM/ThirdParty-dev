if(Threads_FOUND)
  # Save current values
  set(CMAKE_REQUIRED_INCLUDES_TMP ${CMAKE_REQUIRED_INCLUDES})
  set(CMAKE_REQUIRED_LIBRARIES_TMP ${CMAKE_REQUIRED_LIBRARIES})
  # Set flags for building test program
  set(CMAKE_REQUIRED_INCLUDES /usr/include)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

  include(CheckFunctionExists)

  check_function_exists(pthread_setaffinity_np PTHREAD_AFFINITY_LINUX_OK)
  # Restore previous values
  set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES_TMP})
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES_TMP})
endif()
