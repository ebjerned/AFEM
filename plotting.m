clear;
close all;

load results_C12/C12time.txt
load results_C12/C12mutualist.txt
load results_C12/C12prey.txt
load results_C12/C12pred.txt

x = C12mutualist
y = C12prey
z = C12pred
t = C12time

c = 1:numel(t);      %# colors
h = surface([x(:), x(:)], [y(:), y(:)], [z(:), z(:)], ...
    [c(:), c(:)], 'EdgeColor','flat', 'FaceColor','none');
colormap( cool(numel(t)) )
axis tight, grid on, view(55,20)
xlabel("Mutualist")
ylabel("Prey")
zlabel("Predator")
title("Phase diagram of all types")
saveas(gcf, "results_C12/PhaseC12.png")
figure;
load results_sweden/Swedentime.txt
load results_sweden/Swedenmutualist.txt
load results_sweden/Swedenprey.txt
load results_sweden/Swedenpred.txt

x = Swedenmutualist
y = Swedenprey
z = Swedenpred
t = Swedentime

c = 1:numel(t);      %# colors
h = surface([x(:), x(:)], [y(:), y(:)], [z(:), z(:)], ...
    [c(:), c(:)], 'EdgeColor','flat', 'FaceColor','none');
colormap( cool(numel(t)) )
axis tight, grid on, view(55,20)
xlabel("Mutualist")
ylabel("Prey")
zlabel("Predator")
title("Phase diagram of all types")
saveas(gcf, "results_sweden/PhaseSweden.png")


%%
clear;
close all;
plot_solutionPVD("results_C1", "C1solution.pvd", 2, [0, 100, 200, 300, 400, 1000])
plot_solutionPVD("results_C1", "C1solution.pvd", 3, [0, 100, 200, 300, 400, 1000])
plot_solutionPVD("results_C12", "C12solution.pvd", 1, [0, 100, 200, 300, 400, 1000])
plot_solutionPVD("results_C12", "C12solution.pvd", 2, [0, 100, 200, 300, 400, 1000])
plot_solutionPVD("results_sweden", "SwedenSolution.pvd", 1, [0, 50, 100, 600, 1200])
plot_solutionPVD("results_sweden", "SwedenSolution.pvd", 2, [0, 50, 100, 600, 1200])



%% FUNCTIONS FOR PVD-PLOT

function result = plot_solutionPVD(res_folder, pvd_path, component, save_times)
    component_map = ["u", "v", "w"];
    [vtus, time, n] = read_pvd(res_folder + "/" + pvd_path);
    disp(".pvd-file read")
    [p, e, t] = meshVTU(res_folder + "/" + vtus(1));
    disp("Mesh loaded")
    bigdata = assemble_data(res_folder + "/" + vtus);
    disp("Assembled data")
    
    for f = 1:n
        pdeplot(p,e,t, "XYData", bigdata(component,:,f));%, "ZData",  bigdata(2,:,d))
        zlim([0 1])
        if pvd_path == "SwedenSolution.pvd"
            xlim([130 310])
        end
        title(pvd_path + ", Component: " + component_map(component)  + " at time: " + time(f))
        xlabel("x")
        ylabel("y")
        drawnow;
        if sum(ismember(save_times, time(f))) > 0
            saveas(gcf,res_folder + "/" + pvd_path + "_" + component_map(component) + "_"+ time(f) +".png")
        end
    end
    disp("Finished drawing")
end

function [vtus, t, n] = read_pvd(pvd_path)
    dom = xmlread(pvd_path);
    sols = dom.getElementsByTagName("DataSet");
    n = sols.getLength;
    vtus = strings(n, 1);
    t = zeros(n,1);
    for s = 0:n-1
        vtus(s+1) = sols.item(s).getAttributes.getNamedItem("file").getValue;
        t(s+1) = str2num(sols.item(s).getAttributes.getNamedItem("timestep").getValue);
    end
end

function data = assemble_data(path_vec)
    IC = dataVTU(path_vec(1));
    data = zeros(size(IC, 1), size(IC, 2), size(path_vec, 1));
    data(:,:,1) = IC;
    for d = 1:size(path_vec, 1)-1
        data(:,:,d+1) = dataVTU(path_vec(d+1));
    end

end

function [p, e, t] = meshVTU(path)
    dom = xmlread(path);
    p_container = dom.getElementsByTagName("Points").item(0);
    pdata = str2num(p_container.getTextContent)';
    t_container = dom.getElementsByTagName("Cells").item(0).getChildNodes.item(1);
    tdata = str2num(t_container.getTextContent)';
    p = reshape(pdata, 3,length(pdata)/3);
    t = reshape(tdata, 3,length(tdata)/3);
    p = p(1:2,:);
    t = [t; zeros(1, size(t,2))]  + [ones(3, size(t,2)); zeros(1, size(t,2))];
    e = [];
end

function d = dataVTU(path)
    dom = xmlread(path);
    data_container = dom.getElementsByTagName("PointData");
    data = data_container.item(0).getChildNodes.item(1);
    data_size = str2num(data.getAttributes.getNamedItem("NumberOfComponents").getValue);
    pdata = str2num(data.getChildNodes.item(0).getNodeValue)';
    d = reshape(pdata, data_size,length(pdata)/data_size);
end


%% OLD FUNCTIONS
function d = dataXML(path)
    % Renamed
    dom = xmlread(path);
    data_container = dom.getElementsByTagName("PointData");
    data = data_container.item(0).getChildNodes.item(1);
    data_size = str2num(data.getAttributes.getNamedItem("NumberOfComponents").getValue);
    pdata = str2num(data.getChildNodes.item(0).getNodeValue)';
    d = reshape(pdata, data_size,length(pdata)/data_size);
end
function [p, e, t] = meshXML(path)
    % Not used, slower to read mesh in the given xml-format compared to the
    % .vtus
    p = m_meshXML(path, "vertex");
    p = p(2:3, :);
    e = [];
    t = m_meshXML(path, "triangle");
    t = [t(2:4,:); zeros(1, size(t,2))]  + [ones(3, size(t,2)); zeros(1, size(t,2))];
    
    function data_array = m_meshXML(path, data_term)
        dom = xmlread(path);
        data = dom.getElementsByTagName(data_term);

        data_length = data.getLength;
        attrib_length = data.item(0).getAttributes.getLength;
        data_array = zeros(attrib_length, data_length);
        for d = 0:data_length-1
            item = data.item(d);
            attributes = item.getAttributes;
            for a = 0:attrib_length-1
                data_array(a+1, d+1) = str2num(attributes.item(a).getValue);
            end
        end

    end
end




