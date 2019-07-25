function SimulinkExample (block)
% Level-2 MATLAB file S-Function for times two demo.
%   Copyright 1990-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $ 

setup(block);
function setup(block)  
  %% Register number of input and output ports
  block.NumInputPorts  = 1;
  block.NumOutputPorts = 1;

  %% Setup functional port properties to dynamically
  %% inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
  block.InputPort(1).DirectFeedthrough = true;  
  %% Set block sample time to inherited
  block.SampleTimes = [-1 0];
  
  %% Set the block simStateCompliance to default (i.e., same as a built-in block)
  block.SimStateCompliance = 'DefaultSimState';
 
  %% Register methods
  block.RegBlockMethod('Start',@Start);
  block.RegBlockMethod('Outputs',                 @Output);  
  
function Start(block)
 disp('hi')


function Output(block)
% get the input of s-function 
InData = block.InputPort(1).Data;
% get the simulation time:
CurrentTime = get(block,'CurrentTime');

% process the data here here.
OutData = InData * 10;

block.OutputPort(1).Data = OutData;
%%