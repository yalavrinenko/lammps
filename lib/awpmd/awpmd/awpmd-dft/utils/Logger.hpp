//
// Created by cheshire on 06.07.17.
//

#ifndef AWPMD_DFT_LOGGER_HPP
#define AWPMD_DFT_LOGGER_HPP
#include <string>
#include <iostream>

class Logger{
private:
    template <typename ArgT>
    static void m_Out(std::ostream& out, ArgT info_message){
        out << info_message << " ";
    }

    template <class FirstArgT, class ... OtherArgsT>
    static void m_Out(std::ostream& out, FirstArgT info_message, OtherArgsT ... other_args){
        m_Out(out, info_message);
        m_Out(out, other_args...);
    }

public:
    static void SuppressOutput(bool value){
    }

    static void ModuleName(std::string module_name){
    }

    static std::string ModuleName(){
        return "Logger";
    }

    template <class FirstArgT, class ... OtherArgsT>
    static void Info(FirstArgT info_message, OtherArgsT ... other_args){
        std::cout << "["<< "Logger" << " INFO]:";
        m_Out(std::cout, info_message, other_args...);
        std::cout << std::endl;
    }

    template <class FirstArgT, class ... OtherArgsT>
    static void Error(FirstArgT error_message, OtherArgsT ... other_args){
        std::cerr << "\x1B[31m" << "["<< "Logger" << " ERROR]:";
        m_Out(std::cerr, error_message, other_args...);
        std::cerr << "\033[0m" << std::endl;
    }

    template <class FirstArgT, class ... OtherArgsT>
    static void Warning(FirstArgT warning_message, OtherArgsT ... other_args){
        std::cout << "\x1B[34m" << "["<< "Logger" << " WARNING]:";
        m_Out(std::cout, warning_message, other_args...);
        std::cout << "\033[0m" << std::endl;
    }
};

#endif //AWPMD_DFT_LOGGER_HPP
