
# General Todos:

- [x] Create a README.md, as initial documentation.
- [ ] Main PMD features
  - [x] explore ENGINEERING model summary
  - [x] explore ENGINEERING model details 
  - [x] filter elements ENGINEERING model details with a specific rule
  - [x] plot a tree of a DSS file 
  - [x] plot a tree of a PMD ENGINEERING Model 
  - [x] plot a real graph if file has coordinates
  - [ ] explore MATHEMATICAL model summary
  - [x] explore MATHEMATICAL model details
  - [x] filter elements MATHEMATICAL model details with a specific rule
  - [ ] missing elements in MATHEMATICAL model [transformer, switch and shunt]
  - [x] plot a tree of a PMD MATHEMATICAL Model
  - [x] plot a real graph if file has coordinates
  - [ ] explore SIMULATION RESULTS model summary
  - [ ] INVESTIGATE WHY THE LOAD CURRENTS DO NOT MATCH OPENDSS ALTHOUGH VOLTAGES DO -- SEEMS TO BE A BUG
  - [ ] explore SIMULATION RESULTS model details
  - [ ] filter elements SIMULATION RESULTS model details with a specific rule
  - [x] plot a tree of a PMD SIMULATION RESULTS Model with results on top of the graph
  - [ ] if Buchheim assumption fails, then try to split the tree ?!?????????
  - [x] plot a real graph if file has coordinates with results on top of the graph
  - [x] add all test feeders as JLD2


# General Features:

Functional Requirements:

# Basic ENG Model Features:
1. I want to be able to give a PMD ENGINEERING model and get a summary of it and details of its components
2. I want to be able to filter the components of a PMD ENGINEERING model with a specific rule
3. I want to be able to plot a tree of a DSS file
4. I want to be able to plot a tree of a PMD ENGINEERING Model
5. I want to be able to plot a real graph if file has coordinates

# Basic MATH Model Features:
1. I want to be able to give a PMD MATHEMATICAL model and get a summary of it and details of its components
2. I want to be able to filter the components of a PMD MATHEMATICAL model with a specific rule
3. I want to be able to plot a tree of a PMD MATHEMATICAL Model
4. I want to be able to plot a real graph if file has coordinates

# Basic SIM Model Features:
1. I want to be able to give a PMD SIMULATION RESULTS model and get a summary of it and details of its components
2. I want to be able to filter the components of a PMD SIMULATION RESULTS model with a specific rule
3. I want to be able to plot a tree of a PMD SIMULATION RESULTS Model with results on top of the graph
4. I want to be able to plot a real graph if file has coordinates with results on top of the graph

# Advanced PMD Features:
1. I want to be able to compare two models that I have loaded
3. I want to be able edit a model that I have loaded -- specific table item (bus, line, etc.) or add a new item
2. I want to be able to convert a model from one format to another



List the tiles
Tyler.TileProviders.list_providers()