# EGG-Project

This is the source code for the paper titled "An Elite Gene Guided Reproduction Operator for Many-objective Optimization". This java project contains the detailed implementation of EGG to the framework of MOEA/D, NSGA-III, SPEA2-SDE, and $\theta$-DEA.

## Getting Started
These instructions will get you a compy of the project up and running on your local machine for development and testing purposes.

### Prerequisites
#### Java runtime
Java SE Runtime Environment 8 is recommended.

#### jMetal
Our experiment is implemented on the jMetal 4.5 framework, which can be download from (http://jmetal.sourceforge.net/). Using other version of jMetal is not recommended.

### Project Structure and Installation
This project is implemented with Eclipse. So it is recommended that users also use this IDE to run the project.

#### Project structure
There are two important directories in the project.
* [results.] Contain files of the final solutions obtained by each algorithm. When running the experiments, the data will be outputted to a subdirectory named "results/algorithm_name". For example, "results/MOEAD_EGG" saves the results obtained by MOEA/D-EGG algorithm.
* [src.] Contain all source files.

### Installing
1. Build a java project in Eclipse.
2. Unzip this compacted file to the project directory.
3. Refresh the project in Eclipse.
4. Add library "Jama-1.0.2.jar" to the project.

## Running the tests
This section discusses how to reproduce the experiments. The comparison between DE, SBX, EP, and EGG on the framework MOEA/D is illustrated.

1. Run the code to obtain final results of SBX, DE and EGG in different frameworks MOEA/D. The final solutions are saved in the corresponding folder in "results".
2. Obtain results of EP from https://github.com/hxyokokok/EPDE
3. Compare the obtained results by the matlab code "TableS3_MOEAD.m" in the folder 'results'. An example of DTLZ1 problem is provided.

The comparison for the other framework is similar to the above steps.

## Support
If you have any problems, please feel free to connect me via email qlzhu4-c@my.cityu.edu.hk.
