clear all; close all; clc

model = testClass.fromDH([1,2,3,sym('alph1');5,sym('d'),7,8]);

model.T_matrices

model.T(0,2)
