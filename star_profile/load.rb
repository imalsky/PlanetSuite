# load.rb for profiles
require 'scanf'

load '../mesa_dir.rb'

star_test_path = $MESA_DIR + '/star/test'

load star_test_path + '/star_profile/load_star_profile.rb'
