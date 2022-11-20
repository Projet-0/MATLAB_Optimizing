function[s]=fletcher_complet(cx,cy,a,xi,yi,epsilon_fletcher) % d est un doublet cx, cy
    dbtype('fletcher.m') ;

    s_a = [a];
    % On doit d'abord calculer Dk
    while (norm(grad_ctls(cx,cy,xi,yi)) > epsilon_fletcher)

        d = -grad_ctls(cx,cy,xi,yi) ; % Page 13/30
        a = fletcher(cx,cy,a,d,xi,yi) ;

        
        cx = cx + a*d(1) ;
        cy = cy + a*d(2) ;

        s_a = [s_a,a] ;

    end
    
    s = [cx,cy,s_a] ;
end
