#pragma once
#ifndef IO_H
#define IO_H

#include "fields.h"
#include "grid.h"
#include <vector>
#include <string>

void saveCSV(const Fields& F, const std::string& filename);
void saveField3D(const vector<vector<double>>& F, const Grid& G, const string& name);
#endif