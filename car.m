classdef car
    %CAR Class definition of car object.
    %   Object that holds properties of cars and allows for easier
    %   manipulation
    %
    %   carOBJ = CAR(x,y,l,w,v,h) creates a car object with initializations
    %   given by the position (x,y), size (l,w), and pose (v,h).
    
    properties
        X = 0;
        Y = 0;
        velocity = 0;
        heading = 0;
        color = 'yellow'
    end
    
    properties (Access = private)
        L 
        W 
        body
        plots
    end
    
    methods
        %constructor 
        function c = car(x,y,l,w,v,h)
            %CAR Construct an instance of this class
            %   Constructor can take up to 4 arguments to initialize 
            %   position and dimensions of the car.
            if nargin == 0 
                x = 0;
                y = 0;
                l = 1;
                w = 0.5;
                v = 0;
                h = 0;
            end
            import polyshape.*
            c.X         = x;
            c.Y         = y;
            c.velocity  = v;
            c.heading   = h;
            c.L         = l;
            c.W         = w;
            c.body = polyshape([x+l/2 x+l/2 x-l/2 x-l/2],[y+w/2 y-w/2 y-w/2 y+w/2]);
            c.body = rotate(c.body, rad2deg(h), [x y]);
        end
        
        function c = update(c,x,y,v,h)
            %must input full state
            if nargin < 4
                error('4 Inputs required.');
            end
            
            %compute state differences
            xmove = x - c.X;
            ymove = y - c.Y;
            turn  = rad2deg(h - c.heading);
            
            %apply state change
            c.body = translate(c.body, xmove, ymove);
            c.body = rotate(c.body, turn, [x y]);
            
            %update state values
            c.X = x;
            c.Y = y;
            c.velocity = v;
            c.heading  = h;
            %c.color = ;
            
        end
    
       function c = plot(c)
           c.plots = plot(c.body, 'FaceColor', c.color);
       end
       
       %setters
       function c = set.velocity(c, v)
           if v < 0
               warning('Velocity cannot be negative');
           end
           c.velocity = v;
       end
       
       function c = set.heading(c, h)
           if h > 2*pi()
               c.heading = -(h - 2*pi());
           else
               c.heading = h;
           end
           
       end
       
       %getters
       function head = get.heading(c)
           if c.heading == []
               head = 0;
           end
           head = c.heading;
       end
    end
        
end

