#!/bin/bash
make clean
make 
rm -rf png
mkdir  png
./process_data ./params.in
