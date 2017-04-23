# Only run coverage from linux 0.6 build on travis.
get(ENV, "TRAVIS_OS_NAME", "")       == "linux" || exit()
get(ENV, "TRAVIS_JULIA_VERSION", "") == "0.6"   || exit()

using Coverage

Codecov.submit(Codecov.process_folder())
