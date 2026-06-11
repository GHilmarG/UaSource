
function mesh=SelectRefinementEdge(mesh)


d12=(mesh.coordinates(mesh.elements(:,1),1)-mesh.coordinates(mesh.elements(:,2),1)).^2+...
    (mesh.coordinates(mesh.elements(:,1),2)-mesh.coordinates(mesh.elements(:,2),2)).^2;

d23=(mesh.coordinates(mesh.elements(:,2),1)-mesh.coordinates(mesh.elements(:,3),1)).^2+...
    (mesh.coordinates(mesh.elements(:,2),2)-mesh.coordinates(mesh.elements(:,3),2)).^2;

d31=(mesh.coordinates(mesh.elements(:,3),1)-mesh.coordinates(mesh.elements(:,1),1)).^2+...
    (mesh.coordinates(mesh.elements(:,3),2)-mesh.coordinates(mesh.elements(:,1),2)).^2;

[~,K]=max([d12 d23 d31],[],2);

for I=1:size(mesh.elements,1)
    mesh.elements(I,:)= circshift(mesh.elements(I,:),-K(I));
end


end


