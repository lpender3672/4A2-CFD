
#include <iostream>
#include <fstream>
#include <string>
#include "types.h"
#include <cmath>


bool read_settings(const std::string& fpath, t_appvars& av, t_bconds& bcs);

void write_settings(const std::string& fpath, const t_appvars& av, const t_bconds& bcs);
