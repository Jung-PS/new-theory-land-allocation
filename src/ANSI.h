/*
 * ANSI escape sequences for Linux terminal emulators
 * ANSI.h
 * Released under the MIT License, MMC051 Contributor
 */
#pragma once 
#include <iostream> 
#include <string>
namespace ANSI {
      	void escape(){
		std::cout << "\033[";
	}	
	void reset() {
		escape();
		std::cout << "H" << std::flush;
	}
	void color(int id){
		escape();
		std::cout << std::to_string(id) << "m";
	}
	void code(std::string str){
		escape();
		std::cout << str;
	}
}

