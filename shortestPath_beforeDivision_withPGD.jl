using DelimitedFiles
using LinearAlgebra
using Statistics
using Optim


using PyPlot
using PyCall
@pyimport matplotlib.patches as patch

const DEBUG_ORDERWALLS = false
const DEBUG_CENTROID = false
const DEBUG_LENGTH = false
const DEBUG_SHORTESTPATH = false
const DEBUG_EXTRACTVERTICES = false


#const label_n1 = 17
#const label_n2 = 16
#const label_mother = 25
const FACTOR_CENTER = 1.0
const N = 3;
const THRESHOLD_min = 1e-6;
const THRESHOLD_max = 1e10;
const THRESHOLD_angle = 70;

# Read init files
function read_init(source)
    
    # import file and delete comments (and empty lines)
    init_file = readdlm(source, '\t'; header=false, skipstart=0, skipblanks=true, comments=true, comment_char='#')[:,1];

    # Extract number of cells, walls, vertices
    num_cell, num_wall, num_vertex = split(init_file[1]);
    num_cell = parse(Int64, num_cell);
    num_wall = parse(Int64, num_wall);
    num_vertex = parse(Int64, num_vertex);
    delete = popfirst!(init_file)

    # Extract topology matrix
    topo = Array{Int64}(undef,0,5)
    for i = 1:1:num_wall
        n, v1, v2, c1, c2 = split(init_file[1]);
        topo = vcat(topo, [parse(Int64,n) parse(Int64,v1) parse(Int64,v2) parse(Int64,c1) parse(Int64,c2)]);
        delete = popfirst!(init_file); 
    end

    if init_file[1] ==" "
        delete = popfirst!(init_file); 
    end

    # Extract vertices position
    dim = parse(Int64, split(init_file[1])[2]);
    delete = popfirst!(init_file);
    vertices = zeros(num_vertex,dim)
    for i = 1:1:num_vertex
        coord = parse.(Float64, split(init_file[1]));
        vertices[i,:] = coord;
        delete = popfirst!(init_file); 
    end

    if init_file[1] ==" "
        delete = popfirst!(init_file); 
    end

    # Extract wall variables
    dim = parse(Int64, split(init_file[1])[3]);
    delete = popfirst!(init_file);
    walls = zeros(num_wall,dim+1)
    for i = 1:1:num_wall
        vars = parse.(Float64, split(init_file[1]));
        walls[i,:] =  vars;
        delete = popfirst!(init_file); 
    end

    if init_file[1] ==" "
        delete = popfirst!(init_file); 
    end

    # Extract cell variables
    dim = parse(Int64, split(init_file[1])[2]);
    delete = popfirst!(init_file);
    cells = zeros(num_cell,dim)
    for i = 1:1:num_cell
        vars = parse.(Float64, split(init_file[1]));
        cells[i,:] = vars;
        delete = popfirst!(init_file); 
    end

    return topo, vertices, cells

end

# Plot a single cell
function plot_cell(topo, vertices,color, flag, fig,ax)

    for i = 1:1:size(topo,1)
        v1 = topo[i,4];
        v2 = topo[i,5];
        if flag v1 = v1+1 end;
        if flag v2 = v2+1 end;
        x_coord = [vertices[v1,1], vertices[v2,1]];
        y_coord = [vertices[v1,2], vertices[v2,2]];
        #z_coord = [vertices[v1,3], vertices[v2,3]];
        #plot3D(x_coord, y_coord,z_coord,color)
        plot(x_coord, y_coord,color,marker="o")
    end

end


function calculateCentroid(topo, vertices, flag, fig, ax)

   # From: http://coding-experiments.blogspot.com/2009/09/xna-quest-for-centroid-of-polygon.html

   area = 0.0;
   Cx = 0.0;
   Cy = 0.0;

   for i = 1:1:size(topo,1)
    v1 = topo[i,4];
    v2 = topo[i,5];
    x_coord = [vertices[v1,1], vertices[v2,1]];
    y_coord = [vertices[v1,2], vertices[v2,2]];
    DEBUG_CENTROID && print("x = ", x_coord, "\n")
    DEBUG_CENTROID && print("y = ", y_coord, "\n")
    tmp = x_coord[1] * y_coord[2] - x_coord[2] * y_coord[1];
    DEBUG_CENTROID && print("temp = ", tmp, "\n")
    area = area + tmp;
    Cx = Cx + (x_coord[1] + x_coord[2]) *tmp;
    Cy = Cy + (y_coord[1] + y_coord[2]) *tmp;
    DEBUG_CENTROID && print("Cx = ", Cx, "\n")
    DEBUG_CENTROID && print("Cy = ", Cy, "\n")
    flag && ax.plot(x_coord, y_coord,"black", marker = "o")
   end

   area = area/2;
   center_x = Cx/(6*area);
   center_y = Cy/(6*area);
    

    DEBUG_CENTROID && print("center_x = ", center_x, "\n")
    DEBUG_CENTROID && print("center_y = ", center_y, "\n")

    flag && ax.plot(center_x,center_y, marker = "o");

    return center_x, center_y


end


function extract_vertices(topo, vertices, fig, ax)

    new_vertices = Array{Float64}(undef,0,3)

    new_topo = copy(topo);

    unique_vertices = unique(topo[:,4:5]),1;
    map_vertices = fill(-5,size(unique_vertices[1],1));
    map_vertices[:,1] = unique_vertices[1];

    DEBUG_EXTRACTVERTICES && print("map vertices = ", map_vertices, "\n");

    for i = 1:1:size(topo,1)

        map_positionP1 = getindex(findall(x -> x == topo[i,4], map_vertices[:,1]),1);
        map_positionP2 = getindex(findall(x -> x == topo[i,5], map_vertices[:,1]),1);

        new_topo[i,4] = map_positionP1;
        new_topo[i,5] = map_positionP2;

    end

    DEBUG_EXTRACTVERTICES && print("new topo = ", new_topo, "\n");

    new_vertices = Array{Float64}(undef,0,3);
    for i = 1:1:size(unique_vertices[1],1)
        v1 = unique_vertices[1][i] + 1;
        new_vertices = vcat(new_vertices, vertices[v1,1:3]');

    end
    DEBUG_EXTRACTVERTICES && print("new vertices = ", new_vertices, "\n");

    #plot_cell(new_topo, new_vertices,"grey", false, fig,ax)

    return new_topo, new_vertices

end

function best_fitplane(points)

    # Calculate the best fit plane 
    DM = copy(points[:,1:3]);
    DM[:,3] .= 1;
    B = DM\points[:,3]; # z = B[1] x + B[2] y + B[3] equation

    return B

end

function plane_rotation(B, points)

    # Rotate the points to a plane parallel to XY-plane
    cost = 1 / sqrt(B[1]^2 + B[2]^2 + 1);
    sint = sqrt( (B[1]^2 + B[2]^2) / (B[1]^2 + B[2]^2 + 1) );
    u1 = -B[2] / (B[1]^2 + B[2]^2 + 1);
    u2 = B[1] / (B[1]^2 + B[2]^2 + 1);

    r11 = cost+ u1^2*(1-cost);
    r12 = u1*u2* (1-cost);
    r13 = u2 * sint;
    r22 = cost+ u2^2*(1-cost);
    r23 = -u1*sint;
    r33 = cost;

    R = [r11  r12 r13
         r12  r22 r23
        -r13 -r23 r33];

    xy_points = zeros(size(points,1),3);
    for i = 1:1:size(points,1)
        new_point = R * [points[i,1];points[i,2]; points[i,3]];
        
        xy_points[i,1:3] = new_point';
    end

    max = maximum(xy_points, dims = 1);
    min = minimum(xy_points, dims = 1);

    xc = min[1,1] + (max[1,1] - min[1,1])/2;
    yc = min[1,2] + (max[1,2] - min[1,2])/2;

    xy_points[:,1] = xy_points[:,1] .- xc;
    xy_points[:,2] = xy_points[:,2] .- yc;

    return xy_points
end

function extract_new(n1,n2,topo, vertices, fig, ax)
    walls_n1 = findall(x -> x == n1, topo[:,2:3])
    index_n1 = getindex.(walls_n1, [1 2])
    walls_1 = topo[index_n1[:,1],1:5]

    walls_n1 = findall(x -> x== n2,  walls_1[:,2:3])
    index_n1 = getindex.(walls_n1, [1 2])
    new_wall_topo = walls_1[index_n1[:,1],1:5]
    list = [1:1:size(walls_1,1);]
    incomplete_list = filter!(x->!(x in index_n1[:,1]),list)
    walls_end1 = walls_1[incomplete_list,1:5]

    walls_n2 = findall(x -> x == n2, topo[:,2:3])
    index_n2 = getindex.(walls_n2, [1 2])
    walls_2 = topo[index_n2[:,1],1:5]

    walls_n2 = findall(x -> x== n1,  walls_2[:,2:3])
    index_n2 = getindex.(walls_n2, [1 2])
    list = [1:1:size(walls_2,1);]
    incomplete_list = filter!(x->!(x in index_n2[:,1]),list)
    walls_end2 = walls_2[incomplete_list,1:5]

    walls_topo_end = vcat(walls_end1, walls_end2)
    walls_topo_end_ordered = order_walls(walls_topo_end,vertices)

    unique_external = unique(walls_topo_end_ordered[:,4:5]); 
    vertex_external = vertices[unique_external,:]

    unique_new = unique(new_wall_topo[:,4:5]);
    vertex_new = vertices[unique_new,:];

    max_x = maximum(vertex_new[:,1]);
    find_max = getindex.(findall(x->x==max_x, vertices[:,1]),[1])

    find_topo_new = getindex.(findall(x-> x== find_max[1], new_wall_topo[:,4]))
    new_wall_topo_end = vcat(new_wall_topo[find_topo_new,:],new_wall_topo)
    new_wall_export = unique(new_wall_topo_end[:,:], dims = 1);

    #plot_cell(walls_topo_end_ordered, vertices,"black", false, fig,ax)
    #plot_cell(new_wall_topo,vertices,"red", false, fig, ax)


    return new_wall_export, walls_topo_end_ordered, vertex_external

end

function merge_cells(n1,n2,topo, vertices,fig,ax)
    walls_n1 = findall(x -> x == n1, topo[:,2:3])
    index_n1 = getindex.(walls_n1, [1 2])
    walls_1 = topo[index_n1[:,1],1:5]

    list = [1:1:size(walls_1,1);]
    walls_end1 = walls_1[list,1:5]

    walls_n2 = findall(x -> x == n2, topo[:,2:3])
    index_n2 = getindex.(walls_n2, [1 2])
    walls_2 = topo[index_n2[:,1],1:5]

    walls_n2 = findall(x -> x== n1,  walls_2[:,2:3])
    index_n2 = getindex.(walls_n2, [1 2])
    list = [1:1:size(walls_2,1);]
    incomplete_list = filter!(x->!(x in index_n2[:,1]),list)
    walls_end2 = walls_2[incomplete_list,1:5]

    walls_topo_end_disordered = vcat(walls_end1, walls_end2)


    #plot_cell(walls_topo_end_disordered, vertices,"green",true, fig,ax)

    return walls_topo_end_disordered

end


function order_walls(topo,vertices)

    # walls_n = findall(x -> x == n, topo[:,2:3])
    # DEBUG_ORDERWALLS && print("walls_n = ",walls_n, "\n")
    # index = getindex.(walls_n, [1 2])
    # disordered_topo = topo[index[:,1],:];
    disordered_topo =topo
    ordered_topo = Array{Float64}(undef,0,5)
    ordered_topo = vcat(ordered_topo, disordered_topo[1,:]')
    i = 1;
    DEBUG_ORDERWALLS && print("disordered_topo = ",disordered_topo, "\n")
    DEBUG_ORDERWALLS && print("ordered_topo = ",ordered_topo, "\n")

    disordered_topo = disordered_topo[2:end,1:5]
    while !isempty(disordered_topo)
        find_1 = findall(x -> x == ordered_topo[i,5], disordered_topo[:,4])
        find_2 = findall(x -> x == ordered_topo[i,5], disordered_topo[:,5])
        if !isempty(find_1)
            i = i+1;
            index = getindex.(find_1, [1])
            ordered_topo = vcat(ordered_topo, disordered_topo[index[1,1],:]')
            disordered_topo = disordered_topo[1:end .!= index[1,1],1:5]
        elseif !isempty(find_2)
            i = i+1;
            index = getindex.(find_2, [1])
            ordered_topo = vcat(ordered_topo, disordered_topo[index[1,1],:]')
            ordered_topo[i,4] = disordered_topo[index[1,1],5]
            ordered_topo[i,5] = disordered_topo[index[1,1],4]
            disordered_topo = disordered_topo[1:end .!= index[1,1],1:5]
        elseif isempty(find_1) && isempty(find_2)
            #print(ordered_topo);
            print("Warning: Extra elements during merging\n")
            #print(disordered_topo)
            return convert(Array{Int64},ordered_topo)
        end
        #DEBUG_ORDERWALLS && print("disordered_topo = ",disordered_topo, "\n")
        #DEBUG_ORDERWALLS && print("ordered_topo = ",ordered_topo, "\n")
    end
    DEBUG_ORDERWALLS && print("disordered_topo = ",disordered_topo, "\n")
    DEBUG_ORDERWALLS && print("ordered_topo = ",ordered_topo, "\n")

    return convert(Array{Int64},ordered_topo)

end


function shortest_path(label_mother, topo,vertices, cells, topoE, verticesE, cellsE, parents, sourceFib, sourcePGD, sourceLabelsPGD, sourceFibPGD)


    (fig, ax) = subplots(1, 1, figsize=(7,5))

    cell_n3 = findall(x -> x == label_mother, cells[:,2])
    n3 = getindex.(cell_n3, [1])[1] -1

    parent_file = readdlm(parents, ','; header=true)[1];

    find_dauthers = findall(x -> x == label_mother, parent_file[:,2])

    label_n1 = convert(Int64, parent_file[find_dauthers[1],1]);
    label_n2 = convert(Int64,parent_file[find_dauthers[2],1]);

    cell_n1 = findall(x -> x == label_n1, cellsE[:,2])
    n1 = getindex.(cell_n1, [1])[1] -1
    cell_n2 = findall(x -> x == label_n2, cellsE[:,2])
    n2 = getindex.(cell_n2, [1])[1] -1

    DEBUG_SHORTESTPATH && print("label1= ", label_n1, "\n")
    DEBUG_SHORTESTPATH && print("label2= ", label_n2, "\n")
    DEBUG_SHORTESTPATH && print("n1= ", n1, "\n")
    DEBUG_SHORTESTPATH && print("n2= ", n2, "\n")
    DEBUG_SHORTESTPATH && print("n3= ", n3, "\n")

    # Identify the serrounding walls and the points on the new wall from t+dt experiment configuration
    all_topo_merged = merge_cells(n1,n2,topoE, verticesE,fig,ax)

    # extract topology and vertices only of the selected cell at time t+dt
    topo_merged, cell3D_verticesE = extract_vertices(all_topo_merged, verticesE, fig, ax) 

    # coefficient of the plane equation best fit points
    B = best_fitplane(cell3D_verticesE);
    # rotate points in a plane parallel to XY plane
    cell_verticesE = plane_rotation(B, cell3D_verticesE)
    plot_cell(topo_merged, cell_verticesE,"grey", false, fig,ax)

    topo_new_dis, topo_mother, vertices_external = extract_new(n1,n2,topo_merged, cell_verticesE, fig, ax)
    topo_new = order_walls(topo_new_dis,cell_verticesE)

    # Identify the edges of the cell that will divide at time t
    walls_n = findall(x -> x == n3, topo[:,2:3])
    index = getindex.(walls_n, [1 2])
    disordered_topo = topo[index[:,1],:];
    all_topo_cell = order_walls(disordered_topo,vertices)

    # extract topology and vertices only of the selected cell at time t
    topo_cell, cell3D_vertices = extract_vertices(all_topo_cell, vertices, fig, ax) 
    # coefficient of the plane equation best fit points
    B = best_fitplane(cell3D_vertices);
    # rotate points in a plane parallel to XY plane
    cell_vertices = plane_rotation(B, cell3D_vertices)
    #topo_cell = order_walls(topo_cell_dis,cell_vertices)


    # Calculate centroid of cell at time t
    cx, cy = calculateCentroid(topo_cell,cell3D_vertices, false, fig, ax)
    cx = cx*FACTOR_CENTER;
    cy = cy*FACTOR_CENTER;


   # cell_vertices = register_cells(vertices_external[:,1:2], cell_vertices_noreg[:,1:2])
    #plot_cell(topo_cell, cell_vertices,"green", false, fig,ax)
    
    plot_cell(topo_cell, cell_vertices,"grey", false, fig,ax)
    ax.set_aspect("equal")

    # # Calculate fibril orientation
    fibril_vect = compute_angle_v02(label_mother,sourceFib, cx, cy, B, fig, ax, "red")
    m_fibril = (fibril_vect[3,2] - fibril_vect[1,2]) / (fibril_vect[3,1] - fibril_vect[1,1]);


    # Calculate PGD orientation
    #PGD_vect = compute_angle(label_mother,sourcePGD, sourceLabelsPGD, cx, cy, B, fig, ax, "green")
    #m_PGD = (PGD_vect[3,2] - PGD_vect[1,2]) / (PGD_vect[3,1] - PGD_vect[1,1]);
    m_PGD = 10;
 
    # Find tangent to new wall from configuration t+dt
    size_topo_new = round(Int64,size(topo_new,1)/2);
    DEBUG_SHORTESTPATH && print("midd topo new = ", size_topo_new, "\n")
    start_point = size_topo_new-8;
    end_point = size_topo_new+8;
    if start_point<=0 || end_point<=0
        start_point = 1;
        end_point = round(Int64,size(topo_new,1));
    end
    v1 = topo_new[start_point,4]; 
    v2 = topo_new[end_point,4];
    x_coord = [cell_verticesE[v1,1], cell_verticesE[v2,1]];
    y_coord = [cell_verticesE[v1,2], cell_verticesE[v2,2]];
    m_new = (y_coord[2]-y_coord[1]) / (x_coord[2]-x_coord[1]);
    b = y_coord[1] - m_new * x_coord[1];

    # Plot tangent to new wall
    end1 = topo_new[1,4];
    end2 = topo_new[end,4];
    x_coordEnd = [cell_verticesE[end1,1], cell_verticesE[end2,1]];
    y1 = m_new*x_coordEnd[1]+b;
    y2 = m_new*x_coordEnd[2]+b;
    plot([x_coordEnd[1], x_coordEnd[2]], [y1,y2], color = "black", "--")

    # Calculate centroid of cell at time t
    cx, cy = calculateCentroid(topo_cell,cell_vertices, true, fig, ax)
    cx = cx*FACTOR_CENTER;
    cy = cy*FACTOR_CENTER;

    list = [0:1:size(topo_cell,1)-1;]
    DEBUG_SHORTESTPATH && print(list, "\n")

    stored_length = Array{Float64}(undef,0,8)
    min_length = Array{Float64}(undef,0,8)

    #i = 0
    for i = 0:1:size(list,1)-1

        list = [0:1:size(topo_cell,1)-1;]
        if i == 0 
            d = [0, 1, size(topo_cell,1)-1]
        elseif i == size(topo_cell,1)-1
            d = [0, size(topo_cell,1)-2, size(topo_cell,1)-1]
        else
            d = [i-1, i, i+1]
        end

        incomplete_list = filter!(x->!(x in d),list)
        DEBUG_SHORTESTPATH && print("i =", i, "\n")
        DEBUG_SHORTESTPATH && print("incomplete list = ", incomplete_list, "\n")

        v1 = topo_cell[i+1,4];
        v2 = topo_cell[i+1,5];
        x_coord = [cell_vertices[v1,1], cell_vertices[v2,1]];
        y_coord = [cell_vertices[v1,2], cell_vertices[v2,2]];
        m1 = (y_coord[2]-y_coord[1])/ (x_coord[2]-x_coord[1])

        if isinf(m1)
            x_vect = fill(x_coord[1],N);
            step_y = (y_coord[2]-y_coord[1])/(N-1);
            y_vect = [y_coord[1]:step_y:y_coord[2];];
        else
            step_x = (x_coord[2]-x_coord[1])/(N-1);
            x_vect = [x_coord[1]:step_x:x_coord[2];];
            y_vect = m1 .* (x_vect.-x_coord[1]) .+ y_coord[1];
        end

        DEBUG_SHORTESTPATH && plot(x_vect, y_vect,"o", color = "red")
        #j = 2
        for j in incomplete_list

            edge_storage = Array{Float64}(undef,0,8)
            v1 = topo_cell[j+1,4];
            v2 = topo_cell[j+1,5];
            xj_coord = [cell_vertices[v1,1], cell_vertices[v2,1]];
            yj_coord = [cell_vertices[v1,2], cell_vertices[v2,2]];
            m2 = (yj_coord[2]-yj_coord[1])/(xj_coord[2]-xj_coord[1])

            if isinf(m2)
                xj_vect = fill(xj_coord[1],N);
                step_y = (yj_coord[2]-yj_coord[1])/(N-1);
                yj_vect = [yj_coord[1]:step_y:yj_coord[2];];
            else
                step_x = (xj_coord[2]-xj_coord[1])/(N-1);
                xj_vect = [xj_coord[1]:step_x:xj_coord[2];];
                yj_vect = m2 .* (xj_vect.-xj_coord[1]) .+ yj_coord[1];
            end
            
            DEBUG_SHORTESTPATH && plot(xj_vect, yj_vect,"o", color = "blue")

            #z = 3
            for z = 1:1:size(x_vect,1)
               for k = 1:1:size(xj_vect,1)
                    x1 = x_vect[z];
                    y1 = y_vect[z];
                    x3 = xj_vect[k];
                    y3 = yj_vect[k]

                    length,m,theta = calculate_length(x1,y1,cx,cy,x3,y3,m_new,fig,ax,DEBUG_SHORTESTPATH);

                    if !isnan(length)
                     stored_length = vcat(stored_length, [x1 y1 cx cy x3 y3 length theta])
                     edge_storage = vcat(edge_storage, [x1 y1 cx cy x3 y3 length theta])
                    end
               end

            end

            if !isempty(edge_storage)
                edge_storage = sortslices(edge_storage,dims=1, by=x->x[7])
                min_length = vcat(min_length, edge_storage[1,:]')
            end

       end   

    end

    stored_length = sortslices(stored_length,dims=1, by=x->x[7])
    stored_length = round.(stored_length[:,:], digits = 4)

    stored_length = unique(stored_length,dims = 1)

    print("Shortest = ", stored_length[1,7], "\n")
    shortest = stored_length[1,7] *1.02;
    plot_short = true
    k = 1
    short_stored = Array{Float64}(undef,0,8)
    while plot_short==true 
        if (stored_length[k,7] <= shortest)
            short_stored = vcat(short_stored,stored_length[k,:]');
        else
         plot_short = false;
        end
        k = k+1;
    end

    L, m, theta = calculate_length(short_stored[1,1],short_stored[1,2],short_stored[1,3],short_stored[1,4],short_stored[1,5],short_stored[1,6],m_new, fig,ax,true, color_plot = "orange")
    short_stored = sortslices(short_stored,dims=1, by=x->x[8]);
  

    for k = 1:1:size(short_stored,1)
        plot_flag = false;
        if (k ==1) || (k==size(short_stored,1)) 
            plot_flag = true;
        end
        L, m, theta = calculate_length(short_stored[k,1],short_stored[k,2],short_stored[k,3],short_stored[k,4],short_stored[k,5],short_stored[k,6],m_new, fig,ax,plot_flag, color_plot = "blue")
      if plot_flag
        print("k = ", round(k, digits = 4), " L = ", round(L, digits = 4), " theta = ", round(theta, digits = 4), "\n")
      end
    end


    # print("Longest = ", stored_length[end,7], "\n")
    # longest = stored_length[end,7] *0.999;
    # plot_long = true
    # k = size(stored_length,1)
    # long_tangent = Array{Float64}(undef,0,1)
    # while plot_long==true 
    #     if stored_length[k,7] >= longest
    #      L, m, theta = calculate_length(stored_length[k,1],stored_length[k,2],stored_length[k,3],stored_length[k,4],stored_length[k,5],stored_length[k,6],m_new,fig,ax,true, color_plot = "green")
    #      long_tangent = vcat(long_tangent,m);
    #     else
    #         plot_long = false;
    #     end
    #     k = k-1;
    # end
    # theta_longest = rad2deg.(atan.(abs.((m_new.-long_tangent)./(1 .+ long_tangent.*m_new))));

    theta_fibril = rad2deg.(atan.(abs.((m_new.-m_fibril)./(1 .+ m_fibril.*m_new))));
    #theta_PGD = rad2deg.(atan.(abs.((m_new.-m_PGD)./(1 .+ m_PGD.*m_new))));
    #theta_FibPGD = rad2deg.(atan.(abs.((m_fibril.-m_PGD)./(1 .+ m_PGD.*m_fibril))));
    theta_PGD = 100
    theta_FibPGD = 100

    #warning = check_angles(label_mother,theta_FibPGD, sourceFibPGD)
    warning = true

    save = [label_mother stored_length[1,7] stored_length[1,8] short_stored[1,7] short_stored[1,8] theta_fibril theta_PGD theta_FibPGD]


    return stored_length, short_stored, theta_fibril, theta_PGD, theta_FibPGD, warning, save

end



function shortest_pathAllCells(label_mother, topo,vertices, cells, sourceFib, sourceLabelsFib, sourcePGD, sourceLabelsPGD, sourceFibPGD)


    (fig, ax) = subplots(1, 1, figsize=(7,5))

    cell_n3 = findall(x -> x == label_mother, cells[:,2])
    n3 = getindex.(cell_n3, [1])[1] -1

    DEBUG_SHORTESTPATH && print("label1= ", label_n1, "\n")
    DEBUG_SHORTESTPATH && print("n3= ", n3, "\n")

    # Identify the edges of the cell that will divide at time t
    walls_n = findall(x -> x == n3, topo[:,2:3])
    index = getindex.(walls_n, [1 2])
    disordered_topo = topo[index[:,1],:];
    all_topo_cell = order_walls(disordered_topo,vertices)

    # extract topology and vertices only of the selected cell at time t
    topo_cell, cell3D_vertices = extract_vertices(all_topo_cell, vertices, fig, ax) 
    # coefficient of the plane equation best fit points
    B = best_fitplane(cell3D_vertices);
    # rotate points in a plane parallel to XY plane
    cell_vertices = plane_rotation(B, cell3D_vertices)

    # Calculate centroid of cell at time t
    cx, cy = calculateCentroid(topo_cell,cell3D_vertices, false, fig, ax)
    cx = cx*FACTOR_CENTER;
    cy = cy*FACTOR_CENTER;

    plot_cell(topo_cell, cell_vertices,"grey", false, fig,ax)
    ax.set_aspect("equal")

    # # Calculate fibril orientation
    fibril_vect = compute_angle(label_mother,sourceFib, sourceLabelsFib, cx, cy, B, fig, ax, "red")
    m_fibril = (fibril_vect[3,2] - fibril_vect[1,2]) / (fibril_vect[3,1] - fibril_vect[1,1]);

    # Calculate PGD orientation
    PGD_vect = compute_angle(label_mother,sourcePGD, sourceLabelsPGD, cx, cy, B, fig, ax, "green")
    m_PGD = (PGD_vect[3,2] - PGD_vect[1,2]) / (PGD_vect[3,1] - PGD_vect[1,1]);

    # Calculate centroid of cell at time t
    cx, cy = calculateCentroid(topo_cell,cell_vertices, true, fig, ax)
    cx = cx*FACTOR_CENTER;
    cy = cy*FACTOR_CENTER;

    list = [0:1:size(topo_cell,1)-1;]
    DEBUG_SHORTESTPATH && print(list, "\n")

    stored_length = Array{Float64}(undef,0,7)
    min_length = Array{Float64}(undef,0,7)

    for i = 0:1:size(list,1)-1

        list = [0:1:size(topo_cell,1)-1;]
        if i == 0 
            d = [0, 1, size(topo_cell,1)-1]
        elseif i == size(topo_cell,1)-1
            d = [0, size(topo_cell,1)-2, size(topo_cell,1)-1]
        else
            d = [i-1, i, i+1]
        end

        incomplete_list = filter!(x->!(x in d),list)
        DEBUG_SHORTESTPATH && print("i =", i, "\n")
        DEBUG_SHORTESTPATH && print("incomplete list = ", incomplete_list, "\n")

        v1 = topo_cell[i+1,4];
        v2 = topo_cell[i+1,5];
        x_coord = [cell_vertices[v1,1], cell_vertices[v2,1]];
        y_coord = [cell_vertices[v1,2], cell_vertices[v2,2]];
        m1 = (y_coord[2]-y_coord[1])/ (x_coord[2]-x_coord[1])

        if isinf(m1)
            x_vect = fill(x_coord[1],N);
            step_y = (y_coord[2]-y_coord[1])/(N-1);
            y_vect = [y_coord[1]:step_y:y_coord[2];];
        else
            step_x = (x_coord[2]-x_coord[1])/(N-1);
            x_vect = [x_coord[1]:step_x:x_coord[2];];
            y_vect = m1 .* (x_vect.-x_coord[1]) .+ y_coord[1];
        end

        DEBUG_SHORTESTPATH && plot(x_vect, y_vect,"o", color = "red")
        #j = 2
        for j in incomplete_list

            edge_storage = Array{Float64}(undef,0,7)
            v1 = topo_cell[j+1,4];
            v2 = topo_cell[j+1,5];
            xj_coord = [cell_vertices[v1,1], cell_vertices[v2,1]];
            yj_coord = [cell_vertices[v1,2], cell_vertices[v2,2]];
            m2 = (yj_coord[2]-yj_coord[1])/(xj_coord[2]-xj_coord[1])

            if isinf(m2)
                xj_vect = fill(xj_coord[1],N);
                step_y = (yj_coord[2]-yj_coord[1])/(N-1);
                yj_vect = [yj_coord[1]:step_y:yj_coord[2];];
            else
                step_x = (xj_coord[2]-xj_coord[1])/(N-1);
                xj_vect = [xj_coord[1]:step_x:xj_coord[2];];
                yj_vect = m2 .* (xj_vect.-xj_coord[1]) .+ yj_coord[1];
            end
            
            DEBUG_SHORTESTPATH && plot(xj_vect, yj_vect,"o", color = "blue")

            #z = 3
            for z = 1:1:size(x_vect,1)
               for k = 1:1:size(xj_vect,1)
                    x1 = x_vect[z];
                    y1 = y_vect[z];
                    x3 = xj_vect[k];
                    y3 = yj_vect[k]

                    length,m = calculate_length(x1,y1,cx,cy,x3,y3,fig,ax,DEBUG_SHORTESTPATH);

                    if !isnan(length)
                     stored_length = vcat(stored_length, [x1 y1 cx cy x3 y3 length])
                     edge_storage = vcat(edge_storage, [x1 y1 cx cy x3 y3 length])
                    end
               end

            end

            if !isempty(edge_storage)
                edge_storage = sortslices(edge_storage,dims=1, by=x->x[7])
                min_length = vcat(min_length, edge_storage[1,:]')
            end

       end   

    end

    stored_length = sortslices(stored_length,dims=1, by=x->x[7])
    stored_length = round.(stored_length[:,:], digits = 4)

    stored_length = unique(stored_length,dims = 1)

    print("Shortest = ", stored_length[1,7], "\n")
    shortest = stored_length[1,7] *1.02;
    plot_short = true
    k = 1
    short_stored = Array{Float64}(undef,0,7)
    while plot_short==true 
        if (stored_length[k,7] <= shortest)
            short_stored = vcat(short_stored,stored_length[k,:]');
        else
         plot_short = false;
        end
        k = k+1;
    end

    L, m = calculate_length(short_stored[1,1],short_stored[1,2],short_stored[1,3],short_stored[1,4],short_stored[1,5],short_stored[1,6],fig,ax,true, color_plot = "orange")
    short_stored = sortslices(short_stored,dims=1, by=x->x[7]);
  

    for k = 1:1:size(short_stored,1)
        plot_flag = false;
        if (k ==1) || (k==size(short_stored,1)) 
            plot_flag = true;
        end
        L, m = calculate_length(short_stored[k,1],short_stored[k,2],short_stored[k,3],short_stored[k,4],short_stored[k,5],short_stored[k,6], fig,ax,plot_flag, color_plot = "blue")
      if plot_flag
        print("k = ", round(k, digits = 4), " L = ", round(L, digits = 4), "\n")
      end
    end


    theta_fibril = rad2deg.(atan.(abs.((m.-m_fibril)./(1 .+ m_fibril.*m))));
    theta_PGD = rad2deg.(atan.(abs.((m.-m_PGD)./(1 .+ m_PGD.*m))));
    theta_FibPGD = rad2deg.(atan.(abs.((m_fibril.-m_PGD)./(1 .+ m_PGD.*m_fibril))));

    warning = check_angles(label_mother,theta_FibPGD, sourceFibPGD)

    save = [label_mother theta_fibril theta_PGD theta_FibPGD]


    return theta_fibril, theta_PGD, theta_FibPGD, warning, save

end


# Calculate length of circle arc
function calculate_length(x1,y1,x2,y2,x3,y3,fig, ax, flag; color_plot = "black")

    # coordinates mid point first side
    xs = x1 + (x2-x1)/2;
    ys = (y2-y1)/(x2-x1) * (xs-x1) + y1;
    if abs(x1-x2)<=1e-10
        ys = (y2-y1)/2 + y1; 
    end
    DEBUG_LENGTH && print("xs = ", xs, "\n")
    DEBUG_LENGTH && print("ys = ", ys, "\n")
    # coordinates mid point second side
    xss = x2 + (x3-x2)/2;
    yss = (y3-y2)/(x3-x2) * (xss-x2) + y2;
    if abs(x2-x3)<=1e-10
        yss = (y3-y2)/2 + y2; 
    end

    DEBUG_LENGTH && print("xss = ", xss, "\n")
    DEBUG_LENGTH && print("yss = ", yss, "\n")

    # coefficienti rette perpendicolari
    m1 = -(x2-x1)/(y2-y1)
    m2 = -(x3-x2)/(y3-y2)
    DEBUG_LENGTH && print("m1 = ", m1, "\n")
    DEBUG_LENGTH && print("m2 = ", m2, "\n")

    if abs(m1)<THRESHOLD_min
        m1 = 0
    end
    if abs(m2)<THRESHOLD_min
        m2 = 0
    end
    
    DEBUG_LENGTH && ax.plot([x1,x2], [y1,y2], "-",)
    DEBUG_LENGTH && ax.plot([x2,x3], [y2,y3], "-")

    if (isinf(m1) && isinf(m2)) || (iszero(m1) && iszero(m2))
        
        L = sqrt((x2-x1)^2 + (y2-y1)^2) + sqrt((x3-x2)^2 + (y3-y2)^2);
        flag && ax.plot([x1,x3], [y1,y3], "-", color = color_plot)

        m = -(x3-x1)/(y3-y1)
        b = y1-m*x1
        y1 = m*x1+b;
        y3 = m*x3+b;
        flag && ax.plot([x1,x3], [y1 y3], "--")
    
    else

        if isinf(m1)
            x0 = xs;
            y0 = m2 * (x0-xss) + yss;
        elseif isinf(m2)
            x0 = xss;
            y0 = m1 * (x0-xs) + ys;
        else
            # coordinate centro circonferenza
            x0 = (yss-ys-m2*xss + m1* xs) / (m1-m2);
            y0 = m1 * (x0-xs) + ys;
        end

        DEBUG_LENGTH && print("x0 = ", x0, "\n")
        DEBUG_LENGTH && print("y0 = ", y0, "\n")

        # raggio della circonferenza passante per i tre punti
        r = sqrt((x0-x1)^2 + (y0-y1)^2)
        DEBUG_LENGTH && print("r = ", r, "\n")


        if !isinf(r) && (r<THRESHOLD_max)
            
            # Se la somma e' positiva allora i vertici sono in senso orario, altrimenti antiorario
            clockwise = (x2-x1)*(y2+y1) + (x3-x2)*(y3+y2) + (x1-x3)*(y1+y3);
            DEBUG_LENGTH && print("clockwise = ", clockwise, "\n")

            if clockwise > 0
                theta = rad2deg(atan((y3-y0),(x3-x0)))
                alpha = rad2deg(atan((y1-y0),(x1-x0)))
            else
                theta = rad2deg(atan((y1-y0),(x1-x0)))
                alpha = rad2deg(atan((y3-y0),(x3-x0)))
            end

            if theta < 0 
                theta = 360+theta;
            end
            if alpha<0 
                alpha = 360+alpha;
            end

            if alpha >= theta
                beta = alpha-theta;
            elseif theta>alpha
                beta = 360-theta+alpha;
            end

            # angolo compreso arco di circonferenza
            DEBUG_LENGTH && print("alpha = ", alpha, "\n")
            DEBUG_LENGTH && print("theta = ", theta, "\n")
            DEBUG_LENGTH && print("beta = ", beta, "\n")

            m = -(x2-x0)/(y2-y0);
            b = y2-m*x2
            y1 = m*(x2*0.1)+b;
            y3 = m*(x2*10)+b;
            flag && ax.plot([x2*0.1, x2*10], [y1, y3], color = "black", "--", linewidth = 0.1)

           if abs(beta) < THRESHOLD_angle

                # Lunghezza arco di circonferenza = lunghezza cell wall
                L = beta /360 * 2*pi * r;
                DEBUG_LENGTH && print("L = ", L, "\n")


                arcPlot = patch.Arc((x0,y0), 2*r, 2*r, 0, theta, alpha, color = color_plot)
                flag && ax.add_patch(arcPlot)
                DEBUG_LENGTH && ax.plot([xs,x0], [ys,y0], "--")
                DEBUG_LENGTH && ax.plot([xss,x0], [yss,y0], "--")

           else
                L = NaN
           end 

        elseif isinf(r) || (r>THRESHOLD_max)
            L = sqrt((x2-x1)^2 + (y2-y1)^2) + sqrt((x3-x2)^2 + (y3-y2)^2);
            flag && ax.plot([x1, x3], [y1, y3], color = color_plot)

            m = -(x3-x1)/(y3-y1)
            b = y3-m*x3
            y1 = m*x1+b;
            y3 = m*x3+b;
            flag && ax.plot([x1,x3], [y1 y3], color = "black", "--", linewidth = 0.1)
        end

        ax.set_aspect("equal")

    end


    return L, m

end


# Calculate length of circle arc
function calculate_length(x1,y1,x2,y2,x3,y3,m_new,fig, ax, flag; color_plot = "black")

    # coordinates mid point first side
    xs = x1 + (x2-x1)/2;
    ys = (y2-y1)/(x2-x1) * (xs-x1) + y1;
    if abs(x1-x2)<=1e-10
        ys = (y2-y1)/2 + y1; 
    end
    DEBUG_LENGTH && print("xs = ", xs, "\n")
    DEBUG_LENGTH && print("ys = ", ys, "\n")
    # coordinates mid point second side
    xss = x2 + (x3-x2)/2;
    yss = (y3-y2)/(x3-x2) * (xss-x2) + y2;
    if abs(x2-x3)<=1e-10
        yss = (y3-y2)/2 + y2; 
    end

    DEBUG_LENGTH && print("xss = ", xss, "\n")
    DEBUG_LENGTH && print("yss = ", yss, "\n")

    # coefficienti rette perpendicolari
    m1 = -(x2-x1)/(y2-y1)
    m2 = -(x3-x2)/(y3-y2)
    DEBUG_LENGTH && print("m1 = ", m1, "\n")
    DEBUG_LENGTH && print("m2 = ", m2, "\n")

    if abs(m1)<THRESHOLD_min
        m1 = 0
    end
    if abs(m2)<THRESHOLD_min
        m2 = 0
    end
    
    DEBUG_LENGTH && ax.plot([x1,x2], [y1,y2], "-",)
    DEBUG_LENGTH && ax.plot([x2,x3], [y2,y3], "-")

    if (isinf(m1) && isinf(m2)) || (iszero(m1) && iszero(m2))
        
        L = sqrt((x2-x1)^2 + (y2-y1)^2) + sqrt((x3-x2)^2 + (y3-y2)^2);
        flag && ax.plot([x1,x3], [y1,y3], "-", color = color_plot)

        m = -(x3-x1)/(y3-y1)
        b = y1-m*x1
        y1 = m*x1+b;
        y3 = m*x3+b;
        flag && ax.plot([x1,x3], [y1 y3], "--")
    
    else

        if isinf(m1)
            x0 = xs;
            y0 = m2 * (x0-xss) + yss;
        elseif isinf(m2)
            x0 = xss;
            y0 = m1 * (x0-xs) + ys;
        else
            # coordinate centro circonferenza
            x0 = (yss-ys-m2*xss + m1* xs) / (m1-m2);
            y0 = m1 * (x0-xs) + ys;
        end

        DEBUG_LENGTH && print("x0 = ", x0, "\n")
        DEBUG_LENGTH && print("y0 = ", y0, "\n")

        # raggio della circonferenza passante per i tre punti
        r = sqrt((x0-x1)^2 + (y0-y1)^2)
        DEBUG_LENGTH && print("r = ", r, "\n")


        if !isinf(r) && (r<THRESHOLD_max)
            
            # Se la somma e' positiva allora i vertici sono in senso orario, altrimenti antiorario
            clockwise = (x2-x1)*(y2+y1) + (x3-x2)*(y3+y2) + (x1-x3)*(y1+y3);
            DEBUG_LENGTH && print("clockwise = ", clockwise, "\n")

            if clockwise > 0
                theta = rad2deg(atan((y3-y0),(x3-x0)))
                alpha = rad2deg(atan((y1-y0),(x1-x0)))
            else
                theta = rad2deg(atan((y1-y0),(x1-x0)))
                alpha = rad2deg(atan((y3-y0),(x3-x0)))
            end

            if theta < 0 
                theta = 360+theta;
            end
            if alpha<0 
                alpha = 360+alpha;
            end

            if alpha >= theta
                beta = alpha-theta;
            elseif theta>alpha
                beta = 360-theta+alpha;
            end

            # angolo compreso arco di circonferenza
            DEBUG_LENGTH && print("alpha = ", alpha, "\n")
            DEBUG_LENGTH && print("theta = ", theta, "\n")
            DEBUG_LENGTH && print("beta = ", beta, "\n")

            m = -(x2-x0)/(y2-y0);
            b = y2-m*x2
            y1 = m*(x2*0.1)+b;
            y3 = m*(x2*10)+b;
            flag && ax.plot([x2*0.1, x2*10], [y1, y3], color = "black", "--", linewidth = 0.1)

           if abs(beta) < THRESHOLD_angle

                # Lunghezza arco di circonferenza = lunghezza cell wall
                L = beta /360 * 2*pi * r;
                DEBUG_LENGTH && print("L = ", L, "\n")


                arcPlot = patch.Arc((x0,y0), 2*r, 2*r, 0, theta, alpha, color = color_plot)
                flag && ax.add_patch(arcPlot)
                DEBUG_LENGTH && ax.plot([xs,x0], [ys,y0], "--")
                DEBUG_LENGTH && ax.plot([xss,x0], [yss,y0], "--")

           else
                L = NaN
           end 

        elseif isinf(r) || (r>THRESHOLD_max)
            L = sqrt((x2-x1)^2 + (y2-y1)^2) + sqrt((x3-x2)^2 + (y3-y2)^2);
            flag && ax.plot([x1, x3], [y1, y3], color = color_plot)

            m = -(x3-x1)/(y3-y1)
            b = y3-m*x3
            y1 = m*x1+b;
            y3 = m*x3+b;
            flag && ax.plot([x1,x3], [y1 y3], color = "black", "--", linewidth = 0.1)
        end

        ax.set_aspect("equal")

    end

    theta = rad2deg.(atan.(abs.((m_new.-m)./(1 .+ m*m_new))))

    return L, m, theta

end


function compute_angle(n3,source, sourceLabels, cx, cy, B, fig, ax, color)

    angle_file = readdlm(source, ','; header=true)[1];
    findcell = findall(x->(x == n3), angle_file[:,1]);
    alpha = angle_file[findcell[1],2];
    m = tan(deg2rad(alpha))

    label_file = readdlm(sourceLabels, ','; header=true)[1];
    findcell = findall(x->(x == n3), label_file[:,1]);
    if size(findcell,1) == 0
        alpha = 180 - alpha;
    end
        
    
    fib_dir = [ 2*cos(deg2rad(alpha))  2*sin(deg2rad(alpha)) 0
                    0 0 0 
                -2*cos(deg2rad(alpha)) -2*sin(deg2rad(alpha)) 0];


    rotated_dir = plane_rotation(B, fib_dir)
    plot(rotated_dir[:,1], rotated_dir[:,2], color)


    return rotated_dir

end

function compute_angle_v02(n3,source, cx, cy, B, fig, ax, color)

    angle_file = readdlm(source, ','; header=true)[1];
    findcell = findall(x->(x == n3), angle_file[:,1]);
    alpha = angle_file[findcell[1],2];
    m = tan(deg2rad(alpha))
        
    
    fib_dir = [ 2*sin(deg2rad(alpha))  2*cos(deg2rad(alpha)) 0
                    0 0 0 
                -2*sin(deg2rad(alpha)) -2*cos(deg2rad(alpha)) 0];


    rotated_dir = plane_rotation(B, fib_dir)
    plot(rotated_dir[:,1], rotated_dir[:,2], color)


    return rotated_dir

end


function check_angles(n3,theta_FibPGD, source)

    angle_file = readdlm(source, ','; header=true)[1];
    findcell = findall(x->(x == n3), angle_file[:,1]);
    alpha = angle_file[findcell[1],2];
    print("MorphographX Fib-PGD angle =", alpha, "\n")
    warning = false;
    if abs(alpha-theta_FibPGD) > 2
        warning = true
    end

    return warning
end


function register_cells(points1, points2)

    points1_res = fixed_resample(points1[:,1], points1[:,2], 50);
    points2_res = fixed_resample(points2[:,1], points2[:,2], 50);

    R_opt = OPA(points2_res,points1_res,0);

    points2_end = zeros(size(points2,1),2)
    for i = 1:1:size(points2,1)
        points2_end[i,:] = R_opt* points2[i,1:2]
    end

    return points2_end

end

function OPA(X1,X2,centered)

    if centered == 0 
        k1 = size(X1,1);
        k2 = size(X2,1);
        k = minimum([k1,k2]);
        m1 = mean(X1,dims = 1);
        X1c = X1[1:k,:] .- repeat(m1,k);
        m2 = mean(X2,dims = 1);
        X2c = X2[1:k,:] .- repeat(m2,k);
    else
        X1c = copy(X1);
        X2c = copy(X2);
    end

    normX1c=sqrt(tr(transpose(X1c) * X1c));
    normX2c=sqrt(tr(transpose(X2c) * X2c));

    if (normX1c==0 || normX2c==0)
        u=zeros(size(X2c,1),1);
        v=zeros(size(X1c,1),1);
        lambda=0;
    else
        u,lambda,v = svd(transpose(X2c) * X1c /normX1c/normX2c)
    end

    if det(transpose(X2c)*X1c)<0
        u[:,end]=-u[:,end];
        lambda[end,end]=-lambda[end,end];
    end

    R=v*transpose(u);                                
    
    return R

end


function fixed_resample(x::Vector{T}, y::Vector{T}, elperiods::Union{Vector{U}, U}; ) where T<:Float64 where U<:Integer
    

    # initialise resampled arrays as empty
    xᵦ = Vector{Float64}(undef, 0)
    yᵦ = Vector{Float64}(undef, 0)

    for k = 1:1:size(x,1)-1
            x_add  = range(x[k], x[k+1]; length=elperiods+1)[1:end-1]
            y_add =  y[k] .+ (x_add[1:end] .- x[k]).*(y[k+1] .- y[k])./(x[k+1] .- x[k])
            xᵦ = append!(xᵦ, x_add)
            yᵦ = append!(yᵦ, y_add)
    end

    k = 1
    x_add  = range(x[end], x[k]; length=elperiods+1)[1:end-1]
    y_add =  y[k] .+ (x_add[1:end] .- x[k]).*(y[end] .- y[k])./(x[end] .- x[k])
    xᵦ = append!(xᵦ, x_add)
    yᵦ = append!(yᵦ, y_add)


    # last element may be missed, so check and add if needed
    if (xᵦ[end]!= x[end]) && (yᵦ[end]!= y[end])
        append!(xᵦ,x[end])
        append!(yᵦ,y[end])
    end

    return hcat(xᵦ, yᵦ)
end



function mainDividing()

        fileT = ".\\leaf\\analysis\\P5-S5\\P5_S5_24h_parents_12h_MT_output.txt";
        fileTdt = ".\\leaf\\analysis\\P5-S5\\P5_S5_36h_parents_24h_MT_output.txt";
        parents = ".\\leaf\\analysis\\P5-S5\\P5_S5_36h_24h_parents.csv"
        sourceFib = ".\\leaf\\analysis\\P5-S5\\P5_S5_24h_fibrilAngles.csv"
        PGD = ".\\leaf\\analysis\\P5-S5\\P5_S5_24h_AnglesPGD.csv"
        labelsPGD = ".\\leaf\\analysis\\P5-S5\\P5_S5_24h_labelsPGD.csv"
        sourceFibPGD = ".\\leaf\\analysis\\P5-S5\\P5_S5_24h_AnglesDiffFibPGD.csv"
        labels = readdlm(".\\leaf\\analysis\\P5-S5\\P5_S5_24h_36h_labels.csv", ',', header = true)[1]
        outputfile = ".\\leaf\\analysis\\P5-S5\\P5_S5_24h_36h_analysisNew.csv"
        sample_number = 1;
        sample_flag = 1; 

    save_results = Array{Float64}(undef,0,10);

    for z = 1:1:size(labels,1)

        print("--------------------------- \n")
        print("Mother cell label = ", labels[z], "\n")
        topo, vertices, cells = read_init(fileT)
        topoE, verticesE, cellsE = read_init(fileTdt)
        stored_length, short_stored, theta_fibril, theta_PGD, theta_FibPGD, warning, saveLine = shortest_path(labels[z], topo,vertices, cells, topoE, verticesE, cellsE, parents,sourceFib, PGD, labelsPGD, sourceFibPGD)

        print("Fibril angle = ", theta_fibril, "\n")
        print("PGD angle = ", theta_PGD, "\n")
        print("Fibrils-PGD angle = ", theta_FibPGD, "\n")

        if warning == true
            print("Warnings: check angles Fibrils and PGD \n")
        end

        save_results = vcat(save_results, hcat(saveLine, [sample_number sample_flag]))

    end

    writedlm(outputfile,  save_results, ',')


end

function mainAll()

    

    fileT = ".\\leaf\\analysis\\P4-S3\\P4_S3_48h_parents_36h_MT_output.txt";
    sourceFib = ".\\leaf\\analysis\\P4-S3\\P4_S3_48h_fibrilsNew.csv"
    labelsFib =".\\leaf\\analysis\\P4-S3\\P4_S3_48h_labelsFibrils.csv"
    PGD = ".\\leaf\\analysis\\P4-S3\\P4_S3_48h_AnglesPGD.csv"
    labelsPGD = ".\\leaf\\analysis\\P4-S3\\P4_S3_48h_labelsPGD.csv"
    sourceFibPGD = ".\\leaf\\analysis\\P4-S3\\P4_S3_48h_AnglesDiffFibPGD.csv"
    outputfile = ".\\leaf\\analysis\\P4-S3\\P4_S3_48h_analysisAllCells.csv"

    labels = readdlm(".\\leaf\\analysis\\P4-S3\\P4_S3_48h_labelsAll.csv", ',', header = true)[1]

    save_results = Array{Float64}(undef,0,4);


    for z = 1:1:size(labels,1)

        print("--------------------------- \n")
        print("Mother cell label = ", labels[z], "\n")
        topo, vertices, cells = read_init(fileT)
        theta_fibril, theta_PGD, theta_FibPGD, warning, saveLine = shortest_pathAllCells(labels[z], topo,vertices, cells, sourceFib, labelsFib, PGD, labelsPGD, sourceFibPGD)
        print("Fibril angle = ", theta_fibril, "\n")
        print("PGD angle = ", theta_PGD, "\n")
        print("Fibrils-PGD angle = ", theta_FibPGD, "\n")

        if warning == true
            print("Warnings: check angles Fibrils and PGD \n")
        elseif warning == false
            save_results = vcat(save_results, saveLine)
        end

        close("all")

    end

    writedlm(outputfile,  save_results, ',')


end

#mainAll()
mainDividing()

#topo, vertices, cells = read_init(".\\singlecell01.txt")
#shortest_path(topo,vertices,cells)

#(fig, ax) = subplots(1, 1, figsize=(7,5))
# calculate_length(0.66,0.66,1.5,0.5,2.33,0.33,fig,ax,true)
#calculate_length(1.0,0.0,0.5,0.57,0.5,1.5,fig,ax,true)

#--------------------- MARCHANTIA --------------------------------------------------------

# fileT = ".\\analysis\\sample01\\2021_03_17_mSgT_sample01_36h_MT_output.txt";
# fileTdt = ".\\analysis\\sample01\\2021_03_17_mSgT_sample01_48h_output.txt";
# parents = ".\\analysis\\sample01\\sample01_48h_parents.csv"
# sourceFib = ".\\analysis\\sample01\\sample01_36h_fibrilAnglesNew.csv"
# labelsFib =".\\analysis\\sample01\\sample01_36h_labelsFibrils.csv"
# PGD = ".\\analysis\\sample01\\sample01_36h_AnglesPGD.csv"
# labelsPGD = ".\\analysis\\sample01\\sample01_36h_labelsPGD.csv"
# sourceFibPGD = ".\\analysis\\sample01\\sample01_36h_AnglesDiffFibPGD.csv"
# labels = readdlm(".\\analysis\\sample01\\sample01_36h_labels.csv", ',', header = true)[1]
# outputfile = ".\\analysis\\sample01\\sample01_36h_analysisNew.csv"
# sample_number = 1;
# sample_flag = 1; 

# fileT = ".\\analysis\\sample02\\2021_03_17_mSgT_sample02_36h_output.txt";
# fileTdt = ".\\analysis\\sample02\\2021_03_17_mSgT_sample02_48h_output.txt";
# parents = ".\\analysis\\sample02\\sample02_48h_parents.csv"
# sourceFib = ".\\analysis\\sample02\\sample02_36h_fibrilAnglesNew.csv"
# labelsFib =".\\analysis\\sample02\\sample02_36h_labelsFibrils.csv"
# PGD = ".\\analysis\\sample02\\sample02_36h_AnglesPGD.csv"
# labelsPGD = ".\\analysis\\sample02\\sample02_36h_labelsPGD.csv"
# sourceFibPGD = ".\\analysis\\sample02\\sample02_36h_AnglesDiffFibPGD.csv"
# labels = readdlm(".\\analysis\\sample02\\sample02_36h_labels.csv", ',', header = true)[1]
# outputfile = ".\\analysis\\sample02\\sample02_36h_analysisNew.csv"
# sample_number = 1;
# sample_flag = 1; 

# fileT = ".\\analysis\\sample03\\2021_03_17_mSgT_sample03_24h_MT_output.txt";
# fileTdt = ".\\analysis\\sample03\\2021_03_17_mSgT_sample03_36h_output.txt";
# parents = ".\\analysis\\sample03\\sample03_36h_parents.csv"
# sourceFib = ".\\analysis\\sample03\\sample03_24h_fibrilAnglesNew.csv"
# labelsFib =".\\analysis\\sample03\\sample03_24h_labelsFibrils.csv"
# PGD = ".\\analysis\\sample03\\sample03_24h_AnglesPGD.csv"
# labelsPGD = ".\\analysis\\sample03\\sample03_24h_labelsPGD.csv"
# sourceFibPGD = ".\\analysis\\sample03\\sample03_24h_AnglesDiffFibPGD.csv"
# labels = readdlm(".\\analysis\\sample03\\sample03_24h_labels.csv", ',', header = true)[1]
# outputfile = ".\\analysis\\sample03\\sample03_24h_analysisNew.csv"
# sample_number = 1;
# sample_flag = 1; 


# fileT = ".\\analysis\\sample03\\2021_03_17_mSgT_sample03_36h_output.txt";
# fileTdt = ".\\analysis\\sample03\\2021_03_17_mSgT_sample03_48h_output.txt";
# parents = ".\\analysis\\sample03\\sample03_48h_parents.csv"
# sourceFib = ".\\analysis\\sample03\\sample03_36h_fibrilAnglesNew.csv"
# labelsFib =".\\analysis\\sample03\\sample03_36h_labelsFibrils.csv"
# PGD = ".\\analysis\\sample03\\sample03_36h_AnglesPGD.csv"
# labelsPGD = ".\\analysis\\sample03\\sample03_36h_labelsPGD.csv"
# sourceFibPGD = ".\\analysis\\sample03\\sample03_36h_AnglesDiffFibPGD.csv"
# labels = readdlm(".\\analysis\\sample03\\sample03_36h_labels.csv", ',', header = true)[1]
# outputfile = ".\\analysis\\sample03\\sample03_36h_analysisNew.csv"
# sample_number = 1;
# sample_flag = 1; 


# fileT = ".\\analysis\\sample04\\2021_03_17_mSgT_sample04_24h_MT_output.txt";
# fileTdt = ".\\analysis\\sample04\\2021_03_17_mSgT_sample04_36h_MT_output.txt";
# parents = ".\\analysis\\sample04\\sample04_36h_parents.csv"
# sourceFib = ".\\analysis\\sample04\\sample04_24h_fibrilAnglesNew.csv"
# labelsFib =".\\analysis\\sample04\\sample04_24h_labelsFibrils.csv"
# PGD = ".\\analysis\\sample04\\sample04_24h_AnglesPGD.csv"
# labelsPGD = ".\\analysis\\sample04\\sample04_24h_labelsPGD.csv"
# sourceFibPGD = ".\\analysis\\sample04\\sample04_24h_AnglesDiffFibPGD.csv"
# labels = readdlm(".\\analysis\\sample04\\sample04_24h_labels.csv", ',', header = true)[1]
# outputfile = ".\\analysis\\sample04\\sample04_24h_analysisNew.csv"
# sample_number = 1;
# sample_flag = 1; 


# fileT = ".\\analysis\\sample04\\2021_03_17_mSgT_sample04_36h_MT_output.txt";
# fileTdt = ".\\analysis\\sample04\\2021_03_17_mSgT_sample04_48h_output.txt";
# parents = ".\\analysis\\sample04\\sample04_48h_parents.csv"
# sourceFib = ".\\analysis\\sample04\\sample04_36h_fibrilAnglesNew.csv"
# labelsFib =".\\analysis\\sample04\\sample04_36h_labelsFibrils.csv"
# PGD = ".\\analysis\\sample04\\sample04_36h_AnglesPGD.csv"
# labelsPGD = ".\\analysis\\sample04\\sample04_36h_labelsPGD.csv"
# sourceFibPGD = ".\\analysis\\sample04\\sample04_36h_AnglesDiffFibPGD.csv"
# labels = readdlm(".\\analysis\\sample04\\sample04_36h_labels.csv", ',', header = true)[1]
# outputfile = ".\\analysis\\sample04\\sample04_36h_analysisNew.csv"
# sample_number = 1;
# sample_flag = 1; 



# fileT = ".\\analysis\\sample05\\2021_05_27_mSgT_MT_timecourse_sample02_24h_MT_output.txt";
# fileTdt = ".\\analysis\\sample05\\2021_05_27_mSgT_MT_timecourse_sample02_36h_MT_output.txt";
# parents = ".\\analysis\\sample05\\sample05_36h_parents.csv"
# sourceFib = ".\\analysis\\sample05\\sample05_24h_fibrilAnglesNew.csv"
# labelsFib =".\\analysis\\sample05\\sample05_24h_labelsFibrils.csv"
# PGD = ".\\analysis\\sample05\\sample05_24h_AnglesPGD.csv"
# labelsPGD = ".\\analysis\\sample05\\sample05_24h_labelsPGD.csv"
# sourceFibPGD = ".\\analysis\\sample05\\sample05_24h_AnglesDiffFibPGD.csv"
# labels = readdlm(".\\analysis\\sample05\\sample05_24h_labels.csv", ',', header = true)[1]
# outputfile = ".\\analysis\\sample05\\sample05_24h_analysisNew.csv"
# sample_number = 1;
# sample_flag = 1; 


# fileT = ".\\analysis\\sample05\\2021_05_27_mSgT_MT_timecourse_sample02_36h_MT_output.txt";
# fileTdt = ".\\analysis\\sample05\\2021_05_27_mSgT_MT_timecourse_sample02_48h_MT_output.txt";
# parents = ".\\analysis\\sample05\\sample05_48h_parents.csv"
# sourceFib = ".\\analysis\\sample05\\sample05_36h_fibrilAnglesNew.csv"
# labelsFib =".\\analysis\\sample05\\sample05_36h_labelsFibrils.csv"
# PGD = ".\\analysis\\sample05\\sample05_36h_AnglesPGD.csv"
# labelsPGD = ".\\analysis\\sample05\\sample05_36h_labelsPGD.csv"
# sourceFibPGD = ".\\analysis\\sample05\\sample05_36h_AnglesDiffFibPGD.csv"
# labels = readdlm(".\\analysis\\sample05\\sample05_36h_labels.csv", ',', header = true)[1]
# outputfile = ".\\analysis\\sample05\\sample05_36h_analysisNew.csv"
# sample_number = 1;
# sample_flag = 1; 

# fileT = ".\\analysis\\sample05\\2021_05_27_mSgT_MT_timecourse_sample02_48h_MT_output.txt";
# fileTdt = ".\\analysis\\sample05\\2021_05_27_mSgT_MT_timecourse_sample02_60h_MT_output.txt";
# sourceFib = ".\\analysis\\sample05\\sample05_48h_fibrilAngles.csv"
# parents = ".\\analysis\\sample05\\sample05_60h_parents.csv"
# labels = readdlm(".\\analysis\\sample05\\sample05_48h_labels.csv", ',', header = true)[1]
# outputfile = ".\\analysis\\sample05\\sample05_48h_analysis.csv"
# sample_number = 9;
# sample_flag = 1; 

#----------------------------------- LEAF --------------------------------------------

# fileT = ".\\leaf\\analysis\\P5-S5\\P5_S5_0h_MT_output.txt";
# fileTdt = ".\\leaf\\analysis\\P5-S5\\P5_S5_12h_parents_0h_MT_output.txt";
# sourceFib = ".\\leaf\\analysis\\P5-S5\\P5_S5_0h_fibrilAngles.csv"
# parents = ".\\leaf\\analysis\\P5-S5\\P5_S5_12h_0h_parents.csv"
# labels = readdlm(".\\leaf\\analysis\\P5-S5\\P5_S5_0h_labels.csv", ',', header = true)[1]
# outputfile = ".\\leaf\\analysis\\P5-S5\\P5_S5_0h_analysis.csv"
# sample_number = 1;
# sample_flag = 1; 


# fileT = ".\\leaf\\analysis\\P5-S5\\P5_S5_12h_parents_0h_MT_output.txt";
# fileTdt = ".\\leaf\\analysis\\P5-S5\\P5_S5_24h_parents_12h_MT_output.txt";
# sourceFib = ".\\leaf\\analysis\\P5-S5\\P5_S5_12h_fibrilAngles.csv"
# parents = ".\\leaf\\analysis\\P5-S5\\P5_S5_24h_12h_parents.csv"
# labels = readdlm(".\\leaf\\analysis\\P5-S5\\P5_S5_12h_24h_labels.csv", ',', header = true)[1]
# outputfile = ".\\leaf\\analysis\\P5-S5\\P5_S5_12h_24h_analysis.csv"
# sample_number = 2;
# sample_flag = 1; 

# fileT = ".\\leaf\\analysis\\P5-S5\\P5_S5_6h_MT_output.txt";
# fileTdt = ".\\leaf\\analysis\\P5-S5\\P5_S5_24h_parents_6h_MT_output.txt";
# sourceFib = ".\\leaf\\analysis\\P5-S5\\P5_S5_6h_fibrilAngles.csv"
# parents = ".\\leaf\\analysis\\P5-S5\\P5_S5_24h_6h_parents.csv"
# labels = readdlm(".\\leaf\\analysis\\P5-S5\\P5_S5_6h_24h_labels.csv", ',', header = true)[1]
# outputfile = ".\\leaf\\analysis\\P5-S5\\P5_S5_6h_24h_analysis.csv"
# sample_number = 3;
# sample_flag = 1; 

# fileT = ".\\leaf\\analysis\\P5-S5\\P5_S5_24h_parents_12h_MT_output.txt";
# fileTdt = ".\\leaf\\analysis\\P5-S5\\P5_S5_36h_parents_24h_MT_output.txt";
# sourceFib = ".\\leaf\\analysis\\P5-S5\\P5_S5_24h_fibrilAngles.csv"
# parents = ".\\leaf\\analysis\\P5-S5\\P5_S5_36h_24h_parents.csv"
# labels = readdlm(".\\leaf\\analysis\\P5-S5\\P5_S5_24h_36h_labels.csv", ',', header = true)[1]
# outputfile = ".\\leaf\\analysis\\P5-S5\\P5_S5_24h_36h_analysis.csv"
# sample_number = 4;
# sample_flag = 1; 

# fileT = ".\\leaf\\analysis\\P5-S5\\P5_S5_36h_parents_24h_MT_output.txt";
# fileTdt = ".\\leaf\\analysis\\P5-S5\\P5_S5_48h_parents_36h_MT_output.txt";
# sourceFib = ".\\leaf\\analysis\\P5-S5\\P5_S5_36h_fibrilAngles.csv"
# parents = ".\\leaf\\analysis\\P5-S5\\P5_S5_48h_36h_parents.csv"
# labels = readdlm(".\\leaf\\analysis\\P5-S5\\P5_S5_36h_48h_labels.csv", ',', header = true)[1]
# outputfile = ".\\leaf\\analysis\\P5-S5\\P5_S5_36h_48h_analysis.csv"
# sample_number = 5;
# sample_flag = 1; 

# fileT = ".\\leaf\\analysis\\P5-S5\\P5_S5_36h_parents_24h_MT_output.txt";
# fileTdt = ".\\leaf\\analysis\\P5-S5\\P5_S5_48h_parents_36h_MT_output.txt";
# sourceFib = ".\\leaf\\analysis\\P5-S5\\P5_S5_36h_fibrilAngles.csv"
# parents = ".\\leaf\\analysis\\P5-S5\\P5_S5_48h_36h_parents.csv"
# labels = readdlm(".\\leaf\\analysis\\P5-S5\\P5_S5_36h_48h_labels.csv", ',', header = true)[1]
# outputfile = ".\\leaf\\analysis\\P5-S5\\P5_S5_36h_48h_analysis.csv"
# sample_number = 5;
# sample_flag = 1; 

# fileT = ".\\leaf\\analysis\\P5-S5\\P5_S5_48h_parents_36h_MT_output.txt";
# fileTdt = ".\\leaf\\analysis\\P5-S5\\P5_S5_60h_parents_48h_MT_output.txt";
# sourceFib = ".\\leaf\\analysis\\P5-S5\\P5_S5_48h_fibrilAngles.csv"
# parents = ".\\leaf\\analysis\\P5-S5\\P5_S5_60h_48h_parents.csv"
# labels = readdlm(".\\leaf\\analysis\\P5-S5\\P5_S5_48h_60h_labels.csv", ',', header = true)[1]
# outputfile = ".\\leaf\\analysis\\P5-S5\\P5_S5_48h_60h_analysis.csv"
# sample_number = 5;
# sample_flag = 1; 


# fileT = ".\\leaf\\analysis\\P4-S3\\P4_S3_12h_MT_output.txt";
# fileTdt = ".\\leaf\\analysis\\P4-S3\\P4_S3_24h_parents_12h_MT_output.txt";
# parents = ".\\leaf\\analysis\\P4-S3\\P4_S3_24h_12h_parents.csv"
# sourceFib = ".\\leaf\\analysis\\P4-S3\\P4_S3_12h_fibrilsNew.csv"
# labelsFib = ".\\leaf\\analysis\\P4-S3\\P4_S3_12h_labelsFibrils.csv"
# PGD = ".\\leaf\\analysis\\P4-S3\\P4_S3_12h_AnglesPGD.csv"
# labelsPGD = ".\\leaf\\analysis\\P4-S3\\P4_S3_12h_labelsPGD.csv"
# sourceFibPGD = ".\\leaf\\analysis\\P4-S3\\P4_S3_12h_AnglesDiffFibPGD.csv"
# labels = readdlm(".\\leaf\\analysis\\P4-S3\\P4_S3_12h_labels.csv", ',', header = true)[1]
# outputfile = ".\\leaf\\analysis\\P4-S3\\P4_S3_12h_24h_analysisNew.csv"
# sample_number = 1;
# sample_flag = 1; 


# fileT = ".\\leaf\\analysis\\P4-S3\\P4_S3_24h_parents_12h_MT_output.txt";
# fileTdt = ".\\leaf\\analysis\\P4-S3\\P4_S3_36h_parents_24h_MT_output.txt";
# parents = ".\\leaf\\analysis\\P4-S3\\P4_S3_36h_24h_parents.csv"
# sourceFib = ".\\leaf\\analysis\\P4-S3\\P4_S3_24h_fibrilsNew.csv"
# labelsFib = ".\\leaf\\analysis\\P4-S3\\P4_S3_24h_labelsFibrils.csv"
# PGD = ".\\leaf\\analysis\\P4-S3\\P4_S3_24h_AnglesPGD.csv"
# labelsPGD = ".\\leaf\\analysis\\P4-S3\\P4_S3_24h_labelsPGD.csv"
# sourceFibPGD = ".\\leaf\\analysis\\P4-S3\\P4_S3_24h_AnglesDiffFibPGD.csv"
# labels = readdlm(".\\leaf\\analysis\\P4-S3\\P4_S3_24h_labels.csv", ',', header = true)[1]
# outputfile = ".\\leaf\\analysis\\P4-S3\\P4_S3_24h_36h_analysisNew.csv"
# sample_number = 1;
# sample_flag = 1; 

# fileT = ".\\leaf\\analysis\\P4-S3\\P4_S3_36h_parents_24h_MT_output.txt";
# fileTdt = ".\\leaf\\analysis\\P4-S3\\P4_S3_48h_parents_36h_MT_output.txt";
# parents = ".\\leaf\\analysis\\P4-S3\\P4_S3_48h_36h_parents.csv"
# sourceFib = ".\\leaf\\analysis\\P4-S3\\P4_S3_36h_fibrilsNew.csv"
# labelsFib = ".\\leaf\\analysis\\P4-S3\\P4_S3_36h_labelsFibrils.csv"
# PGD = ".\\leaf\\analysis\\P4-S3\\P4_S3_36h_AnglesPGD.csv"
# labelsPGD = ".\\leaf\\analysis\\P4-S3\\P4_S3_36h_labelsPGD.csv"
# sourceFibPGD = ".\\leaf\\analysis\\P4-S3\\P4_S3_36h_AnglesDiffFibPGD.csv"
# labels = readdlm(".\\leaf\\analysis\\P4-S3\\P4_S3_36h_labels.csv", ',', header = true)[1]
# outputfile = ".\\leaf\\analysis\\P4-S3\\P4_S3_36h_48h_analysisNew.csv"
# sample_number = 1;
# sample_flag = 1; 

# fileT = ".\\leaf\\analysis\\P4-S3\\P4_S3_48h_parents_36h_MT_output.txt";
# fileTdt = ".\\leaf\\analysis\\P4-S3\\P4_S3_60h_parents_48h_MT_output.txt";
# sourceFib = ".\\leaf\\analysis\\P4-S3\\P4_S3_48h_fibrils.csv"
# parents = ".\\leaf\\analysis\\P4-S3\\P4_S3_60h_48h_parents.csv"
# labels = readdlm(".\\leaf\\analysis\\P4-S3\\P4_S3_48h_labels.csv", ',', header = true)[1]
# outputfile = ".\\leaf\\analysis\\P4-S3\\P4_S3_48h_60h_analysis.csv"
# sample_number = 5;
# sample_flag = 1; 


# fileT = ".\\leaf\\analysis\\P1-S6\\P1_S6_12h_output.txt";
# fileTdt = ".\\leaf\\analysis\\P1-S6\\P1_S6_24h_output.txt";
# parents = ".\\leaf\\analysis\\P1-S6\\P1_S6_24h_12h_parents.csv"
# sourceFib = ".\\leaf\\analysis\\P1-S6\\P1_S6_12h_fibrilsNew.csv"
# labelsFib = ".\\leaf\\analysis\\P1-S6\\P1_S6_12h_labelsFibrils.csv"
# PGD = ".\\leaf\\analysis\\P1-S6\\P1_S6_12h_AnglesPGD.csv"
# labelsPGD = ".\\leaf\\analysis\\P1-S6\\P1_S6_12h_labelsPGD.csv"
# sourceFibPGD = ".\\leaf\\analysis\\P1-S6\\P1_S6_12h_AnglesDiffFibPGD.csv"
# labels = readdlm(".\\leaf\\analysis\\P1-S6\\P1_S6_12h_labels.csv", ',', header = true)[1]
# outputfile = ".\\leaf\\analysis\\P1-S6\\P1_S6_12h_24h_analysisNew.csv"
# sample_number = 1;
# sample_flag = 1; 

# fileT = ".\\leaf\\analysis\\P1-S6\\P1_S6_24h_output.txt";
# fileTdt = ".\\leaf\\analysis\\P1-S6\\P1_S6_36h_output.txt";
# parents = ".\\leaf\\analysis\\P1-S6\\P1_S6_36h_24h_parents.csv"
# sourceFib = ".\\leaf\\analysis\\P1-S6\\P1_S6_24h_fibrilsNew.csv"
# labelsFib = ".\\leaf\\analysis\\P1-S6\\P1_S6_24h_labelsFibrils.csv"
# PGD = ".\\leaf\\analysis\\P1-S6\\P1_S6_24h_AnglesPGD.csv"
# labelsPGD = ".\\leaf\\analysis\\P1-S6\\P1_S6_24h_labelsPGD.csv"
# sourceFibPGD = ".\\leaf\\analysis\\P1-S6\\P1_S6_24h_AnglesDiffFibPGD.csv"
# labels = readdlm(".\\leaf\\analysis\\P1-S6\\P1_S6_24h_labels.csv", ',', header = true)[1]
# outputfile = ".\\leaf\\analysis\\P1-S6\\P1_S6_24h_36h_analysisNew.csv"
# sample_number = 1;
# sample_flag = 1; 


# fileT = ".\\leaf\\analysis\\P1-S6\\P1_S6_36h_output.txt";
# fileTdt = ".\\leaf\\analysis\\P1-S6\\P1_S6_48h_output.txt";
# sourceFib = ".\\leaf\\analysis\\P1-S6\\P1_S6_36h_fibrils.csv"
# parents = ".\\leaf\\analysis\\P1-S6\\P1_S6_48h_36h_parents.csv"
# labels = readdlm(".\\leaf\\analysis\\P1-S6\\P1_S6_36h_labels.csv", ',', header = true)[1]
# outputfile = ".\\leaf\\analysis\\P1-S6\\P1_S6_36h_48h_analysis.csv"
# sample_number = 5;
# sample_flag = 1; 

# fileT = ".\\leaf\\analysis\\P1-S6\\P1_S6_48h_output.txt";
# fileTdt = ".\\leaf\\analysis\\P1-S6\\P1_S6_60h_output.txt";
# parents = ".\\leaf\\analysis\\P1-S6\\P1_S6_60h_48h_parents.csv"
# sourceFib = ".\\leaf\\analysis\\P1-S6\\P1_S6_48h_fibrilsNew.csv"
# labelsFib = ".\\leaf\\analysis\\P1-S6\\P1_S6_48h_labelsFibrils.csv"
# PGD = ".\\leaf\\analysis\\P1-S6\\P1_S6_48h_AnglesPGD.csv"
# labelsPGD = ".\\leaf\\analysis\\P1-S6\\P1_S6_48h_labelsPGD.csv"
# sourceFibPGD = ".\\leaf\\analysis\\P1-S6\\P1_S6_48h_AnglesDiffFibPGD.csv"
# labels = readdlm(".\\leaf\\analysis\\P1-S6\\P1_S6_48h_labels.csv", ',', header = true)[1]
# outputfile = ".\\leaf\\analysis\\P1-S6\\P1_S6_48h_60h_analysisNew.csv"
# sample_number = 1;
# sample_flag = 1; 