function opi = get_opi(Latitude,scan_node_type);
%============================== Form Orbit Phase ==========================

opi = NaN(length(Latitude),1);

% Descending
% From equator to S. Pole
kd = (scan_node_type == 1 & Latitude <= 0);
[n,e,bin] = histcounts(Latitude(kd),[-90:2:0]);
opi(kd) = abs(bin-46);
% From N. Pole to equator
kd = (scan_node_type == 1 & Latitude >= 0);
[n,e,bin] = histcounts(Latitude(kd),[0:2:90]);
opi(kd) = abs(bin-46)+135;

% Ascending (no index jumps, so only one histcounts)
ka = (scan_node_type == 0);
[n,e,bin] = histcounts(Latitude(ka),[-90:2:90]);
opi(ka) = bin + 45;


