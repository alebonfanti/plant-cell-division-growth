#using CSV
using DelimitedFiles
using Revise
#using DataFrames
using PyPlot
"""
 From the file generated in Morphographx extract data points and connections
     n -> extract the number of point (centers + junctions)
     data_points -> numbering, x, y, z, -1 (it belongs to more than one cell)/lable from segmentation -> Save points and junctions in array
     points_c -> numbering, x, y, z, lable from segmentation -> Save the centers of cells in array
     raw_connections -> connections of each point
"""

function uploadfile(filedir)

    n = parse(Int64, readline(filedir));
    temp_file = string(@__DIR__, "/temp_file.txt");
    raw_points = hcat(zeros(n,5),fill("e",n,1));
    count_c = 0;
    open(filedir) do f
          io_cn = open(temp_file, "w")
           i = 1;
           for l in eachline(f)
              if i <= n+1
                if i != 1
                  a, b, c, d, e, f = split(l);
                  new_raw = [parse(Int64,a) parse(Float64,b) parse(Float64,c) parse(Float64,d) parse(Float64,e) String(f)];
                  raw_points[i-1,:] = new_raw[:];
                  if f=="c" 
                    count_c+=1 
                  end
                end
                i += 1;
              else
                println(io_cn, l);
              end
           end
          close(io_cn);
     end
    raw_connections = readdlm(temp_file);
    points_c = copy(raw_points[1:count_c,1:5])
    rm(temp_file)
     
   
    return n, count_c, raw_points, raw_connections, points_c
end

"""
 Function to plot the geometry
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
Function to remove intermediate segments. The cell wall are aproximated as streight
"""
function streightwall(raw_connections,n_c)

    new_connection = fill(-1,n_c,size(raw_connections)[2])
    new_connection[1:end,1] = raw_connections[1:n_c,1]

    for i = 1:1:n_c
        j = 3
        new = []
        for j=3:1:raw_connections[i,2]+2
            n_point = raw_connections[raw_connections[i,j]+1,2]
            if n_point == 6
              new = push!(new,raw_connections[i,j])
            end
        end
        new_connection[i,2] = size(new,1)
        new_connection[i,3:size(new)[1]+2] = new[1:end]

    end

    max_colum = maximum(new_connection[:,2])+2
    output_connections = copy(new_connection[:,1:max_colum])

    return output_connections
end


"""
 Create init file for tissue
 Functions that generates the topology of each cell: edge_number, cell1, cell2, point1, point2
 each line of edge_cell_connection represents a cell wall (edge_number)
 cell1 and cell2 are the two cells that the edge divides
 point1 and point2 are the two verteces of the wall
"""
function init_cell(n_c, connections, all_connections)
  edge_cell_connection = [0 0 0 0 0]

  for i = 2:1:n_c
    count_edge = 0
    cell2 = 0
    # check if there is only one wall in the cell - then it skips the cell (degenerated cell)
    if (connections[i,2] != 2) || (connections[i,2] != 1)
      for j = 3:1:connections[i,2]+2
          cell1 = connections[i,1]-1
          if j == (connections[i,2]+2)
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
                # In the connections of point p1 look for cell centers (their identifier must be less than number of points)
                index_centers = findall(x->x<n_c, all_connections[p1+1,3:(all_connections[p1+1,2]+2)])
                # Identify the centers of the cell in p1 connection array
                cell_centers = all_connections[p1+1,3:(all_connections[p1+1,2]+2)][index_centers]

                z = 1
                cell2_found = false
                while (cell2_found == false && (z<=size(cell_centers,1)))
                    if (!isempty(findall(x->x==p2, all_connections[cell_centers[z]+1,3:(all_connections[cell_centers[z]+1,2]+2)])) && cell_centers[z] != connections[i,1])
                      cell2 = cell_centers[z]-1
                      cell2_found = true
                    end
                    z = z+1
                end

                edge_cell_connection = vcat(edge_cell_connection,[count_edge cell1 cell2 p1 p2])
                count_edge = count_edge +1
              end
      end
        #print("i ", i, "/", size(connections,1), "\n")
    end
  end
    
    # Remove first line - it is there to avoid errors for the first cell
    edge_cell_connection_return = fill(0,size(edge_cell_connection,1)-1,5)
    edge_cell_connection_return[1:end,:] = edge_cell_connection[2:end,:]
    # Enumerate edges from 0 to their size
    for i=1:size(edge_cell_connection_return,1) edge_cell_connection_return[i,1]=(i-1) end


    return edge_cell_connection_return

end

"""
Calculate the length of the edges
"""
function calculate_length(edge_cell_connection,data_points)
    # Calculate the length of each cell wall
    wall_length = fill(0.0,size(edge_cell_connection,1),2)
    for i = 1:1:size(edge_cell_connection,1)
      p1 = edge_cell_connection[i,4]
      p2 = edge_cell_connection[i,5]
      coord1 = data_points[p1+1,2:4]
      coord2 = data_points[p2+1,2:4]
      wall_length[i,2] = ((coord1[1]-coord2[1])^2 + (coord1[2]-coord2[2])^2 + (coord1[3]-coord2[3])^2)^(0.5)
    end
    return wall_length
end


"""
Re-enumerate cells and verteces to have sequential values
"""

function renumerate_cell_vertex(edge_cell_connection,data_points)

      # Enumarte cells from 0 and substitute the new numbers into the topology matrix
      unique_cells = unique(edge_cell_connection[:,2:3])
      #print(unique_cells)
      cells = fill(0.0,size(unique_cells,1),3)
      cells[:,1] = range(-1, size(unique_cells,1)-2; step = 1)
      cells[:,2] = sort(unique_cells)
      # keep the correspondent lable in morphographx
      for i = 1:1:size(unique_cells,1)
        cells[i,3] = data_points[convert(Int64,cells[i,2])+2,5];
      end
      
      edge_connection_new = zeros(size(edge_cell_connection,1), 2)
      for i = 1:1:size(edge_cell_connection,1)
        C = findall(x->x==edge_cell_connection[i,2], cells[:,2])
        edge_connection_new[i,1] = cells[C[1],1]
        C = findall(x->x==edge_cell_connection[i,3], cells[:,2])
        edge_connection_new[i,2] = cells[C[1],1]
      end

      edge_cell_connection[:,2:3] = edge_connection_new[:,1:2];

      # if edge_cell_connection[i,3]>edge_cell_connection[end,2]
      #   edge_cell_connection[i,3] = -1
      # end

      # Enumarte verteces from 0 and substitute the new numbers into the topology matrix
      unique_vertex = unique(edge_cell_connection[:,4:5])
      vertex = fill(0.0,size(unique_vertex,1),5)
      vertex[:,1] = range(0, size(unique_vertex,1)-1; step = 1)
      vertex[:,2] = unique_vertex
      for i = 1:1:size(unique_vertex,1)
        vertex[i,3:5] = data_points[convert(Int64,vertex[i,2])+1,2:4];
      end

      for i = 1:1:size(edge_cell_connection,1)
        C = findall(x->x==edge_cell_connection[i,4], vertex[:,2])
        edge_cell_connection[i,4] = vertex[C[1],1]
        C = findall(x->x==edge_cell_connection[i,5], vertex[:,2])
        edge_cell_connection[i,5] = vertex[C[1],1]
      end

      return edge_cell_connection, vertex, cells
end


function celldivision(cells, num_var_cells, directoryPath)

  cell_data = fill(0.0,size(cells,1),num_var_cells)
  cell_data[:,36] = cells[:,3]; # Morphographx lables


  if directoryPath != ""

    array = readdlm(directoryPath, ',', header = true)
    dataraw = array[1]

      for i = 1:1:size(dataraw,1)
        
        if dataraw[i,2] != 1

          pos = findall(x->x==dataraw[i,1], cells[:,3])
          if !isempty(pos)
            cell_data[pos[1],37] = dataraw[i,2]-1;
            cell_data[pos[1],38] = 1;
          end

        end
        
      end
  end

  return cell_data

end

function cellselected(cells, num_var_cells, directoryPath)

  cell_data = fill(0.0,size(cells,1),num_var_cells)
  cell_data[:,2] = cells[:,3]; # Morphographx lables


  if directoryPath != ""

    array = readdlm(directoryPath, ',', header = true)
    dataraw = array[1]

      for i = 1:1:size(dataraw,1)
        

          pos = findall(x->x==dataraw[i,1], cells[:,3])
          if !isempty(pos)
            cell_data[pos[1],3] = dataraw[i,1]-1;
            cell_data[pos[1],4] = 1;
          end

  
        
      end
  end

  return cell_data

end


# Write input file for tissue

function print_file(directoryPath_tissue, edge_cell_connection,vertex, cells, wall_length, num_var_cells, num_var_edge, vertex_dim, cell_data)


    wall_data = fill(0.0,size(wall_length,1),num_var_edge)
    wall_data[:,1] = wall_length[:,2];

    io_tissue = open(directoryPath_tissue, "w")
    println(io_tissue, "# Input file generated with Alessandra's script")

    println(io_tissue, size(cell_data,1)-1, " ", size(edge_cell_connection,1), " ", size(vertex,1), "# NumCell NumWall NumVertex")
    for i=1:1:size(edge_cell_connection,1)
        writedlm(io_tissue, edge_cell_connection[i,:]', ' ')
    end
    println(io_tissue, " ")

    println(io_tissue, size(vertex,1), " ", vertex_dim, " ", "# NumVertex Dimension")
    for i=1:1:size(vertex,1)
        writedlm(io_tissue, round.(vertex[i,3:5]',digits=2), ' ')
    end

    println(io_tissue, " ")

    println(io_tissue, size(wall_data,1), " ", "1 ", num_var_edge-1, "# NumEdge NumLength NumVariable")
    for i=1:1:size(wall_data,1)
        writedlm(io_tissue, round.(wall_data[i,:]',digits=2), ' ')
    end

    println(io_tissue, " ")

    println(io_tissue, size(cell_data,1)-1," ", num_var_cells, "# Num_cells Num_variable_cells")
    for i=2:1:size(cell_data,1)
        writedlm(io_tissue, cell_data[i,:]', ' ')
    end

    close(io_tissue)
end


#-------------------------------------------------------------------------------
function main(name)
      print("@ Working directory: ")
      
      println(@__DIR__)

      print("@ File name: ", name, "\n")
      directoryPath = name #string(@__DIR__, name)
      n, n_c, data_points, raw_connections, points_c = uploadfile(directoryPath);
      print("@ Points coordinates and connections uploaded \n")

      #new_connections = streightwall(raw_connections,n_c)
      new_connections = raw_connections;
      print("@ Streight connections created \n")

      plotgeometry(new_connections,data_points,n_c)
      #plotgeometry(raw_connections,raw_points)
      #print("Geometry plotted \n")
      #print(n_c)
      edge_cell_connection = init_cell(n_c, new_connections, raw_connections);
      print("@ Edges matrix created \n")

      wall_length = calculate_length(edge_cell_connection,data_points)
      print("@ Edge length calculated \n")

      edge_cell_connection,vertex,cells = renumerate_cell_vertex(edge_cell_connection,data_points)
      print("@ Cells and verteces re-enumerated \n")


      num_var_cells = 4
      num_var_edge = 3
      vertex_dim = 3
      name_tissue = string(name,"_output.txt") #string(@__DIR__,name[1:end-4],"_output.txt")

      #directoryProliferation = string(@__DIR__, "\\2021_03_17_mSgT_MT_timecourse\\sample01_36h_MT_parents.csv");
      directoryProliferation = ""
      cell_data = cellselected(cells, num_var_cells, directoryProliferation)
      print("@ Cells data matrix ready \n")

      print_file(name_tissue, edge_cell_connection,vertex, cells, wall_length, num_var_cells, num_var_edge, vertex_dim, cell_data)
      print("@ Init file ready! :) \n")

end

#-------------------------------------------------------------------------------
name = "G:\\Laptop_SLCU\\Documents\\Morpho_Tissue_tools\\Morphographx_to_Tissue\\2021_03_17_mSgT_sample01_24h_patch02.txt"
main(name)