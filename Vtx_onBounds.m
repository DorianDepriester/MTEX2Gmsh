function [ vtx_IDs ] = Vtx_onBounds(gB)
    outerBoundary_id=any(gB.grainId==0,2);      % IDs of the faces neighbouring no other grains
    list_IDs=gB.F(outerBoundary_id(:),:);       % Corresponding vertices
    vtx_IDs=unique(list_IDs);
end