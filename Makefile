CXX=g++                                                                       
                                                                                
OBDIR=/home/heinen/anaconda3/envs/openbabel                                   
                                                                                
INCLUDE=-I/home/heinen/anaconda3/envs/openbabel/include/openbabel-2.0/       
LIBS=-L/home/heinen/anaconda3/envs/openbabel/lib                             
                                                                                
CXX_FLAGS=-std=c++11 -O3 -march=native -Wall -D_GLIBCXX_USE_CXX11_ABI=0 
LINKER_FLAGS=-lopenbabel                                                      
                                                                                
all: mcAmon                                                                     
                                                                                
mcAmon: src/mcAmon.cpp src/utils.hpp                                    
	$(CXX) $(INCLUDE) $(LIBS) $(CXX_FLAGS) src/mcAmon.cpp -o mcAmon $(LINKER_FLAGS)
                                                                                
clean:                                                                          
	rm -f mcAmon
	rm *.xyz
