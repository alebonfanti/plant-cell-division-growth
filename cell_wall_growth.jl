using CSV
using DelimitedFiles
using DataFrames
using PyPlot
using LsqFit

"""
 From the file generated in Morphographx extract data points and connections
     n -> extract the number of point (centers + junctions)
     data_points -> numbering, x, y, z, -1 (it belongs to more than one cell)/lable from segmentation -> Save points and junctions in array
     points_c -> numbering, x, y, z, lable from segmentation -> Save the centers of cells in array
     raw_connections -> connections of each point
"""

function uploadfile(filedir)

    n = parse(Int64, readline(filedir));
    print(n)
    temp_file = "G://UCL-Experiments//force_measurement//2022_04_05//camera//m1//1st_stretch//temp_file.txt"
    raw_points = hcat(zeros(n,5));
    count_c = 0;

    open(filedir) do f
          io_cn = open(temp_file, "w")
           i = 1;
           for l in eachline(f)
              if i <= n+1
                if i != 1
                  a, b, c, d, e, f = split(l);
                  new_raw = [parse(Int64,a) parse(Float64,b) parse(Float64,c) parse(Float64,d) parse(Float64,e)];
                  raw_points[i-1,:] = new_raw[:];
                  if f=="c" count_c+=1 end
                end
                i += 1;
              else
                println(io_cn, l);
              end
           end
          close(io_cn);
     end
    raw_connections = readdlm(temp_file);
    centre_connections = raw_connections[1:count_c,:];
    max_col = maximum(raw_connections[count_c+1:end,2])
    edge_connections = raw_connections[count_c+1:end,1:max_col]
    points_c = copy(raw_points[1:count_c,1:5])
    rm(temp_file)

    return n, count_c, raw_points, centre_connections, edge_connections, raw_connections, points_c
end

"""
 Function to plot the geometry directly for the imported file
"""
function plotgeometry(connections, coordinates, n_c)
    fig = figure()
    for i = 1:1:n_c
      if connections[i,2] != 0
          for j = 1:1: connections[i,2]-1
            n1_point = connections[i,j+2]
            n2_point = connections[i,j+3]
            x_coords = [coordinates[n1_point+1,2], coordinates[n2_point+1,2]]
            y_coords = [coordinates[n1_point+1,3], coordinates[n2_point+1,3]]
            z_coords = [coordinates[n1_point+1,4], coordinates[n2_point+1,4]]
            plot3D(x_coords, y_coords,z_coords,"black")
          end
          n1_point = connections[i,3]
          n2_point = connections[i,connections[i,2]+2]
          x_coords = [coordinates[n1_point+1,2], coordinates[n2_point+1,2]]
          y_coords = [coordinates[n1_point+1,3], coordinates[n2_point+1,3]]
          z_coords = [coordinates[n1_point+1,4], coordinates[n2_point+1,4]]
          plot3D(x_coords, y_coords,z_coords,"black")
     end
    end

end


"""
Function to create the cell wall geometry: n -> cell wall number
                                           L1 -> lable of the first cell
                                           L2 -> lable of the second cell (-1 for the edge cells)
                                           P1 -> first point of the edge
                                           P2 -> second point of the edge



 Create init file for tissue
 Functions that generates the topology of each cell: edge_number, cell1, cell2, point1, point2
 each line of edge_cell_connection represents a cell wall (edge_number)
 cell1 and cell2 are the two cells that the edge divides
 point1 and point2 are the two verteces of the wall
"""
function init_cellWall(n_c, connections, all_connections, data_points)
 
    edge_cell_connection = zeros(1,5)
    count_edge = 0
  for i = 1:1:size(connections,1)
    cell2 = 0
    for j = 3:1:connections[i,2]+2
        cell1 = connections[i,1]
        if j == (connections[i,2]+2) # last element of the connection array must be linked to the first one
          p1 = connections[i,connections[i,2]+2]
          p2 = connections[i,3]
        else
          p1 = connections[i,j]
          p2 = connections[i,j+1]
        end

        b1 = [p1;p2]
        b2 = [p2;p1]
        sum_find1 = sum([edge_cell_connection[f,4:5] == b1 for f=1:size(edge_cell_connection,1)])
        sum_find2 = sum([edge_cell_connection[f,4:5] == b2 for f=1:size(edge_cell_connection,1)])

          if (sum_find1+sum_find2) == 0
              # In the connections of point p1 look for cell centers (their identifier must be less than number of center points)
              index_centers = findall(x->x<n_c, all_connections[p1+1,3:(all_connections[p1+1,2]+2)])
              # Identify the centers of the cell in p1 connection array
              cell_centers = all_connections[p1+1,3:(all_connections[p1+1,2]+2)][index_centers]

              z = 1
              cell2_found = false
              while (cell2_found == false && (z<=size(cell_centers,1)))
                  if (!isempty(findall(x->x==p2, all_connections[cell_centers[z]+1,3:(all_connections[cell_centers[z]+1,2]+2)])) && cell_centers[z] != connections[i,1])
                    cell2 = cell_centers[z]
                    cell2_found = true
                  end
                  z = z+1
              end

              if cell2_found == false
                cell2 = -1
              end

              edge_cell_connection = vcat(edge_cell_connection,[count_edge cell1 cell2 p1 p2])
              count_edge = count_edge +1
            end
      end

    end
    edge_cell_connection_return = fill(0,size(edge_cell_connection,1)-1,5)
    edge_cell_connection_return[1:end,:] = edge_cell_connection[2:end,:]

    # Substitute cell number with the Morphographx lable


    for i = 1:1:size(edge_cell_connection_return,1)
        for j = 2:1:3
            if edge_cell_connection_return[i,j] != -1
                edge_cell_connection_return[i,j] = data_points[edge_cell_connection_return[i,j]+1,5]
            end
        end
    end


    return edge_cell_connection_return

end
   
"""
Function to plot the cell walls from the connectivity cell wall geometry (faster)
"""

function plotcellWalls(edge_cell_connection, coordinates, color)
    for i = 1:1:size(edge_cell_connection,1)
            n1_point = edge_cell_connection[i,4]
            n2_point = edge_cell_connection[i,5]
            x_coords = [coordinates[n1_point+1,2], coordinates[n2_point+1,2]]
            y_coords = [coordinates[n1_point+1,3], coordinates[n2_point+1,3]]
            z_coords = [coordinates[n1_point+1,4], coordinates[n2_point+1,4]]
            plot3D(x_coords, y_coords,z_coords,color)
    end

end

"""
Calculate cell wall length
"""
function cellWall_length(edge_cell_connection, data_points)

    unique_cells = unique(edge_cell_connection[:,2:3], dims = 1);
    
    length_walls = zeros(size(unique_cells,1), 3);
    length_walls[:,1:2] = unique_cells[:,1:2];

    for i = 1:1:size(unique_cells,1)
        tot_length = 0;
        for j = 1:1:size(edge_cell_connection,1)
          if ((edge_cell_connection[j,2]==unique_cells[i,1]) && (edge_cell_connection[j,3]==unique_cells[i,2]))  
            # print("i = ", i, "j = ", j, "\n")
            x1 = data_points[edge_cell_connection[j,4]+1,2]
            y1 = data_points[edge_cell_connection[j,4]+1,3]
            z1 = data_points[edge_cell_connection[j,4]+1,4]

            x2 = data_points[edge_cell_connection[j,5]+1,2]
            y2 = data_points[edge_cell_connection[j,5]+1,3]
            z2 = data_points[edge_cell_connection[j,5]+1,4]

            update_length = ((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)^(0.5)

            tot_length = tot_length + update_length
          end
        end
        length_walls[i,3] = tot_length
    end

    return length_walls
end

function loadingFile(name)

    n, n_c, data_points, centre_connections, edge_connections, raw_connections, points_c = uploadfile(directoryPath);
    print("Points coordinates and connections uploaded \n")

    edge_cell_connection = init_cellWall(n_c, centre_connections, raw_connections, data_points)
    print("Cell wall connectivity matrix calculated \n")

    length_walls = cellWall_length(edge_cell_connection,data_points)

    return data_points, edge_cell_connection, length_walls


end


function calculateStrain(walls_old, walls_new, parents) 
  
  walls_new_parents = copy(walls_new)

    for i =1:1:size(walls_new,1)
        if walls_new[i,1] != -1
            parent_pos = findfirst(x -> x == walls_new[i,1], parents[:,1])
            if parent_pos[1] != []
             walls_new_parents[i,1] = parents[parent_pos[1],2]
            end
        end
        if walls_new[i,2] != -1
            parent_pos = findfirst(x -> x == walls_new[i,2], parents[:,1])
            walls_new_parents[i,2] = parents[parent_pos[1],2]
        end
    end
    
    strain = zeros(size(walls_old,1),4);
    strain[:,1:2] = walls_old[:,1:2];


    for i =1:1:size(walls_old,1)
        for j = 1:1:size(walls_new_parents,1)
            if (walls_old[i,1]==walls_new_parents[j,1] && walls_old[i,2]==walls_new_parents[j,2]) || (walls_old[i,1]==walls_new_parents[j,2] && walls_old[i,2]==walls_new_parents[j,1])
                strain[i,3] = strain[i,3]+ walls_new_parents[j,3]
            end
        end
    end

    strain[:,4] = (strain[:,3] .- walls_old[:,3]) ./ walls_old[:,3] * 100; # 1: label 1 old; 2: label 2 old; 3: length; 4: strain

    return strain

end


function plotstrain(edge_cell_connection, coordinates, strains, min_val, max_val)

    (fig, ax) = subplots(1, 1, figsize=(10,10))
    ax.set_facecolor((0, 0, 0))

    min_strain = minimum(strain[:,4]);
    print("Min strain= ", min_strain, "\n")
    max_strain = maximum(strain[:,4]);
    print("Max strain = ", max_strain)

    for i = 1:1:size(edge_cell_connection,1)

        n1_point = edge_cell_connection[i,4]
        n2_point = edge_cell_connection[i,5]
        x_coords = [coordinates[n1_point+1,2], coordinates[n2_point+1,2]]
        y_coords = [coordinates[n1_point+1,3], coordinates[n2_point+1,3]]
        z_coords = [coordinates[n1_point+1,4], coordinates[n2_point+1,4]]

        j = 1;
        found = false
        strain_plot = 0
        index_plot = 0
        while found==false
            if (edge_cell_connection[i,2]==strain[j,1] && edge_cell_connection[i,3]==strain[j,2]) || ((edge_cell_connection[i,3]==strain[j,1]) && (edge_cell_connection[i,2]==strain[j,2]))
                strain_plot = strain[j,4];
                index_plot = round(2/(max_val-min_val) * strain_plot - (max_val + min_val) / (max_val-min_val), digits = 3)
                found = true;
            end
            j = j+1;
        end

        

        if (index_plot < -1) 
          r = 0
          g = 0
          b = 1
        elseif (index_plot >= -1) && (index_plot < -0.5)
          r = 0 
          g = 2*index_plot + 2
          b = 1
        elseif (index_plot >= -0.5) && (index_plot < 0)
          r = 0
          g = 1
          b = -2* index_plot
        elseif (index_plot >= 0) && (index_plot < 0.5)
          r = 2 * index_plot
          g = 1
          b = 0
        elseif (index_plot >= 0.5) && (index_plot < 1)
          r = 1
          g = -2 * index_plot + 2
          b = 0
        elseif (index_plot >= 1)
          r = 1
          g = 0
          b = 0
        end

        if (edge_cell_connection[i,2]==1) || (edge_cell_connection[i,3]==1) 
          r= 0.5  
          g = 0.5
          b = 0.5
        end


         ax.plot(x_coords, y_coords,color = (r,g,b), linewidth=3)
    end
    ax.set(aspect="equal")


end


function wallYMStrain(wallYM, strain, parents, namesave)

  cellwalls = fill(NaN, (size(strain,1),5));
  cellwalls[:,1:4] = strain[:,:];

  wallYM_parents = copy(wallYM);
  for i =1:1:size(wallYM_parents,1)
        #print("label1 = ", wallYM[i,1], " label2 = ", wallYM[i,2], "\n")
        parent_pos = findfirst(x -> x == wallYM[i,1], parents[:,1])
        wallYM_parents[i,1] = parents[parent_pos[1],2]
        parent_pos = findfirst(x -> x == wallYM[i,2], parents[:,1])
        wallYM_parents[i,2] = parents[parent_pos[1],2]
  end

  for i =1:1:size(cellwalls,1)
   for j = 1:1:size(wallYM_parents,1)
      if (cellwalls[i,1]==wallYM_parents[j,1] && cellwalls[i,2]==wallYM_parents[j,2]) || (cellwalls[i,1]==wallYM_parents[j,2] && cellwalls[i,2]==wallYM_parents[j,1])
        if !isnan(cellwalls[i,5])  
          cellwalls[i,5] = (cellwalls[i,5] + wallYM_parents[j,3])/2;
        else
          cellwalls[i,5] = wallYM_parents[j,3];
        end
      end
    end
  end
  #show(stdout, "text/plain",cellwalls)
  cellwalls_noNaN = Array{Float64}(undef,0,5)

  for i = 1:1:size(cellwalls,1)
    if !isnan(cellwalls[i,5]) && (cellwalls[i,3] > -10)
      cellwalls_noNaN = vcat(cellwalls_noNaN, cellwalls[i,:]')
    end
  end

  @. model(x, p) = p[1] *x + p[2]

  p0 = [1e6, 1e6]
  x = vec(cellwalls_noNaN[:,4]);
  y = vec(cellwalls_noNaN[:,5]);
  fit = curve_fit(model, x, y, p0)

  params = fit.param
  print("\n", "Fitting parameters = ", params, "\n")
  xmin = minimum(x);
  xmax = maximum(x);

  xvect = range(xmin, stop=xmax, length=100);
  yvect = model(xvect, params);

  (fig, ax) = subplots(1, 1, figsize=(7,5))
  #fig.suptitle("Sample02", fontsize = 11)
  ax.plot(cellwalls[:,4],cellwalls[:,5], "o")
  ax.plot(xvect,yvect, linewidth = 3)
  ax.grid("on", alpha=0.4, which="major")
  ax.set_xlabel(L"Strain, $\epsilon$", fontsize=11,fontname = "serif")
  ax.set_ylabel("Young's modulus, E (Pa)", fontsize=11,fontname = "serif")
  #ax.set_ylim([0.0, 1.0])

  (fig2, ax2) = subplots(1, 1, figsize=(7,5))
  ax2.plot(cellwalls[:,4],cellwalls[:,3], "o")
  ax2.grid("on", alpha=0.4, which="major")
  ax2.set_ylabel("Length", fontsize=11,fontname = "serif")
  ax2.set_xlabel(L"Strain, $\epsilon$", fontsize=11,fontname = "serif")

  # cellwalls array -> 1: label 1 old; 2: label 2 old; 3: length; 4: strain; 5: YM

  writedlm( namesave,  cellwalls, ',')

  return cellwalls

end

function wallYMRatio(YM1, cellwalls)

  ratio = fill(NaN, (size(YM1,1),4))
  ratio[:,1:2] = YM1[:,1:2];


  for i =1:1:size(ratio,1)
    for j = 1:1:size(cellwalls,1)
       if (ratio[i,1]==cellwalls[j,1] && ratio[i,2]==cellwalls[j,2]) || (ratio[i,1]==cellwalls[j,2] && ratio[i,2]==cellwalls[j,1])
           ratio[i,3] = cellwalls[j,5] / YM1[i,3] ;
           ratio[i,4] = cellwalls[j,4] #strain
       end
     end
   end
 
  # show(stdout, "text/plain",ratio)

  (fig, ax) = subplots(1, 1, figsize=(7,5))
  ax.plot(ratio[:,4],ratio[:,3], "o")
  ax.grid("on", alpha=0.4, which="major")
  ax.set_xlabel(L"Strain, $\epsilon$", fontsize=11,fontname = "serif")
  ax.set_ylabel("Ratio YM", fontsize=11,fontname = "serif")
  

  # average = filter(!isnan, ratio[:,3])
  # (fig2, ax2) = subplots(1, 1, figsize=(7,5))
  # ax2.boxplot(average, # Each column/cell is one box
  #       notch=true, # Notched center
  #       whis=0.95, # Whisker length as a percent of inner quartile range
  #       widths=0.05, # Width of boxes
  #       sym="rs") # Symbol color and shape (rs = red square)

  return ratio

end


#-------------MARCHANTIA 1 PLASMALISATION---------------------------------------------
# sample01; sample02; sample03

# # Mesh segmentation first time point
# name = "\\2021_03_17_mSgT_dev_exp08\\segmentation\\2021_03_17_mSgT_sample01_24h_area01.txt"
# println(@__DIR__)
# directoryPath = string(@__DIR__, name)

# data_points_24, edge_cell_connections_24, length_walls_24 = loadingFile(name);

# # Mesh segmentation second time point
# name = "\\2021_03_17_mSgT_dev_exp08\\segmentation\\2021_03_17_mSgT_sample01_48h_area01.txt"
# directoryPath = string(@__DIR__, name)

# data_points_48,  edge_cell_connections_48, length_walls_48 = loadingFile(name);

# # Parents first-second time point
# name = "\\2021_03_17_mSgT_dev_exp08\\segmentation\\2021_03_17_mSgT_sample01_parents.csv"
# directoryPath = string(@__DIR__, name)
# parents = readdlm(directoryPath, ',', header = true)[1];

# # Strain calculation
# strain = calculateStrain(length_walls_24, length_walls_48, parents);
# plotstrain(edge_cell_connections_24, data_points_24, strain, 40,130)

# # Calculate strain - YM correlation
# name = "\\2021_03_17_mSgT_dev_exp08\\AFM\\2021_03_17_mSgT_sample01_cellwall_YM_line.csv"
# directoryPath = string(@__DIR__, name)
# wallYM = readdlm(directoryPath, ',', header = false);

# namefile = "\\2021_03_17_mSgT_dev_exp08\\2021_03_17_mSgT_sample01_analysis_area01.csv"
# directoryPath = string(@__DIR__, namefile)
# cellwalls = wallYMStrain(wallYM, strain, parents, directoryPath);


#-------------LEAF---------------------------------------------

# change: sample01; sample02, sample05

# Mesh segmentation first time point
# name = "G://UCL-Experiments//force_measurement//2022_04_05//camera//m1//1st_stretch//0_area01_p2.txt"
# println(@__DIR__)
# directoryPath = name #string(@__DIR__, name)

# data_points_24, edge_cell_connections_24, length_walls_24 = loadingFile(name);

# # Mesh segmentation second time point
# name = "G://UCL-Experiments//force_measurement//2022_04_05//camera//m1//1st_stretch//10_area01_p2.txt"
# directoryPath = name # string(@__DIR__, name)

# data_points_48,  edge_cell_connections_48, length_walls_48 = loadingFile(name);

# # Parents first-second time point
# name = "G://UCL-Experiments//force_measurement//2022_04_05//camera//m1//1st_stretch//parent_0_10.csv"
# directoryPath = name #string(@__DIR__, name)
# parents = readdlm(directoryPath, ',', header = true)[1];

# # Strain calculation
# strain = calculateStrain(length_walls_24, length_walls_48, parents);
# plotstrain(edge_cell_connections_24, data_points_24, strain, -15,30)

# Calculate strain - YM correlation
# name = "\\2021_03_04_pmTomato_leaf_7days\\AFM\\2021_03_04_leaf_sample05_F_500_7days_cellwall_YM_line.csv"
# directoryPath = string(@__DIR__, name)
# wallYM = readdlm(directoryPath, ',', header = false);

# namefile = "\\2021_03_04_pmTomato_leaf_7days\\2021_03_04_leaf_sample05_7days_analysis.csv"
# directoryPath = string(@__DIR__, namefile)
# cellwalls = wallYMStrain(wallYM, strain, parents, directoryPath);

# # show(stdout, "text/plain",cellwalls)

# # Calculate correlation YM first-second time point
# name = "\\sample01\\2020_11_05_sample01_area01_F_300_24h_YM_line.csv"
# directoryPath = string(@__DIR__, name)
# wallYM1 = readdlm(directoryPath, ',', header = false);

# YMratio = wallYMRatio(wallYM1, cellwalls);

# show(stdout, "text/plain",YMratio)
# show(stdout,"text/plain", YMratio[findall(x -> x>2.0, YMratio[:,3]),:])
# show(stdout,"text/plain", YMratio[findall(x -> x<0.85, YMratio[:,3]),:])



#-------------MARCHANTIA (double plasmalisation AFM) ------------------------------

# change: sample01; sample04

# #Mesh segmentation first time point
# name = "\\2020_11_05_mSgT_dev_exp07\\segmentation\\sample04\\2020_11_05_mSgT_AFM_sample04_24h_3um.txt"
# println(@__DIR__)
# directoryPath = string(@__DIR__, name)

# data_points_24, edge_cell_connections_24, length_walls_24 = loadingFile(name);

# # Mesh segmentation second time point
# name = "\\2020_11_05_mSgT_dev_exp07\\segmentation\\sample04\\2020_11_05_mSgT_AFM_sample04_48h_3um.txt"
# directoryPath = string(@__DIR__, name)

# data_points_48,  edge_cell_connections_48, length_walls_48 = loadingFile(name);

# # Parents first-second time point
# name = "\\2020_11_05_mSgT_dev_exp07\\segmentation\\sample04\\2020_11_05_sample04_24_48_parents.csv"
# directoryPath = string(@__DIR__, name)
# parents = readdlm(directoryPath, ',', header = true)[1];

# # Strain calculation
# strain = calculateStrain(length_walls_24, length_walls_48, parents);
# plotstrain(edge_cell_connections_24, data_points_24, strain, -10,120)

# # Calculate strain - YM correlation
# name = "\\2020_11_05_mSgT_dev_exp07\\AFM\\Analysis\\2020_11_06_sample04_area01_F_300_48h_YM_line.csv"
# directoryPath = string(@__DIR__, name)
# wallYM = readdlm(directoryPath, ',', header = false);

# namefile = "\\2020_11_05_mSgT_dev_exp07\\2020_11_05_mSgT_sample04_48h_analysis.csv"
# directoryPath = string(@__DIR__, namefile)
# cellwalls = wallYMStrain(wallYM, strain, parents, directoryPath);

# # show(stdout, "text/plain",cellwalls)

# # Calculate correlation YM first-second time point
# name = "\\sample01\\2020_11_05_sample01_area01_F_300_24h_YM_line.csv"
# directoryPath = string(@__DIR__, name)
# wallYM1 = readdlm(directoryPath, ',', header = false);

# YMratio = wallYMRatio(wallYM1, cellwalls);

# show(stdout, "text/plain",YMratio)
# show(stdout,"text/plain", YMratio[findall(x -> x>2.0, YMratio[:,3]),:])
# show(stdout,"text/plain", YMratio[findall(x -> x<0.85, YMratio[:,3]),:])


#----------------------------------------------------------------------



# Mesh segmentation first time point
name = "G:\\Laptop_SLCU\\Documents\\Marchantia_Experiments\\AFM_Confocal_Timecourse\\fig\\2021_03_17_mSgT_sample01_24h_patch02.txt"
println(@__DIR__)
directoryPath = name #string(@__DIR__, name)

data_points_24, edge_cell_connections_24, length_walls_24 = loadingFile(name);

# Mesh segmentation second time point
name = "G:\\Laptop_SLCU\\Documents\\Marchantia_Experiments\\AFM_Confocal_Timecourse\\fig\\2021_03_17_mSgT_sample01_48h_patch02.txt"
directoryPath = name # string(@__DIR__, name)

data_points_48,  edge_cell_connections_48, length_walls_48 = loadingFile(name);

# Parents first-second time point
name = "G:\\Laptop_SLCU\\Documents\\Marchantia_Experiments\\AFM_Confocal_Timecourse\\fig\\sample01_area_left_parents.csv"
directoryPath = name #string(@__DIR__, name)
parents = readdlm(directoryPath, ',', header = true)[1];

# Strain calculation
strain = calculateStrain(length_walls_24, length_walls_48, parents);
plotstrain(edge_cell_connections_24, data_points_24, strain, 40,115)

# Calculate strain - YM correlation
# name = "\\2021_03_04_pmTomato_leaf_7days\\AFM\\2021_03_04_leaf_sample05_F_500_7days_cellwall_YM_line.csv"
# directoryPath = string(@__DIR__, name)
# wallYM = readdlm(directoryPath, ',', header = false);

# namefile = "\\2021_03_04_pmTomato_leaf_7days\\2021_03_04_leaf_sample05_7days_analysis.csv"
# directoryPath = string(@__DIR__, namefile)
# cellwalls = wallYMStrain(wallYM, strain, parents, directoryPath);

# # show(stdout, "text/plain",cellwalls)

# # Calculate correlation YM first-second time point
# name = "\\sample01\\2020_11_05_sample01_area01_F_300_24h_YM_line.csv"
# directoryPath = string(@__DIR__, name)
# wallYM1 = readdlm(directoryPath, ',', header = false);

# YMratio = wallYMRatio(wallYM1, cellwalls);

# show(stdout, "text/plain",YMratio)
# show(stdout,"text/plain", YMratio[findall(x -> x>2.0, YMratio[:,3]),:])
# show(stdout,"text/plain", YMratio[findall(x -> x<0.85, YMratio[:,3]),:])

print("\n", "End")
