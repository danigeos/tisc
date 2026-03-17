# tAo Frequently Asked Questions

### Q1. How can I learn to use this code on my own?
**A1.** After download, follow the instructions in the `README.md` file in the root directory in order to compile. Then go to `tao/demo` and try the examples there to ensure that the code is working properly. Then read carefully the documentation and `tao/doc/template.PRM` and try to experiment modifying the `*.PRM` files of the demos. Then learn how to build `*.UNIT` input files from the examples. Make sure you did all this before contacting the author!

### Q2. I have compilation errors, is there a problem in the program?
**A2.** I tested tAo only in some UNIX operative systems such as Linux (Ubuntu), macOS, and Solaris. Most probably you have a problem with the compilation options or with the libraries. Try to ask for help among your colleagues to succeed in the compilation.

### Q3. I compiled the program, but I never managed to run it, and I always get 'segmentation fault' errors.
**A3.** This means you have some incompatibility between some of the functions called in tAo and the ones available in the libraries of your operative system. You have to edit the code and search for those problems. You will need the help of someone around you who has experience in programming and/or in compiling C programs.

### Q4. Can you help me to solve the problems I have to run tAo?
**A4.** In general, I try to answer the questions you send by email, but usually I cannot say much more than what you just read in this file. Unfortunately, I cannot provide you a close technical support on the compilation process, but this has a good side: you get this program for free. Because my interests are in the scientific use of tAo and its application to geological scenarios, I will be pleased to help you with any problem you have during the development of a model and will be pleased to establish a scientific collaboration.

### Q5. Is it possible to stop a fault and then reactivate it again?
**A5.** You can stop a fault at a given time by using the `time_stop` parameter in the `*.UNIT` file. To reactivate it, repeat the same fault geometry in another `*.UNIT` file with the velocity you wish. Use `tao -hu` to get the possible parameters in every `*.UNIT` file.
