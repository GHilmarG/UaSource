classdef PointerHandle < handle
    properties (SetAccess = private)
        pointer;
        name;
        owner;
    end
    methods
        function obj = PointerHandle(qptr, name, owner)
            obj.owner   = 0;
            obj.pointer = 0;
            obj.name    = 'unknown';
            if nargin>=1
                if isa(qptr, 'PointerHandle')
                    obj.pointer = qptr.pointer;
                    obj.name    = qptr.name;
                    return;
                else
                    obj.pointer = qptr;
                end
                if nargin>=2 & isstr(name)
                    obj.name    = name;
                end
                if nargin>=3
                    obj.owner   = owner;
                end
            end
        end
        function delete(obj)
            if obj.owner
                pointer_handle(obj, 'delete');
            end
        end
    end
end
