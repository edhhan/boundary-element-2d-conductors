% function [G0,H0]=calc_GH0(p0,elem)
%
% Calcul analytique de G0 et H0
%
% En entr�e:
%
% p0   : les coordonn�es x et x dy point d'observation
% elem : la liste des �l�ments sur lesquels porte le calcul
%
% En sortie:
%
% G0 et H0 : vecteurs lignes contenant une valeur d'int�grale par �l�ment

function [G0v,H0v]=calc_GH0(p0,elem)

% Cr�ation des vecteurs de sortie
N=length(elem);
G0v=zeros(1,N);
H0v=zeros(1,N);
seuil=5*eps;

% Boucle sur tous les �l�ments de la structure de donn�es
for i=1:N
	% Donn�es importantes, r�cup�r�es de la structure de donn�es
	p1=elem(i).p1;
	p2=elem(i).p2;
	n =elem(i).n;
	ln=elem(i).ln;

	% Ceci doit �tre calcul� pour chaque point d'observation
	r1=p1-p0;		% Vecteurs r1 et r2, dans le rep�re original
	r2=p2-p0;
	D=dot(r1,n);	% Distance perpendiculaire � l'�l�ment (peut �tre n�gative, et c'est correct!)
	d=D*n;			% Vecteur distance (orient� avec la normale � l'�l�ment)
	L1=dot(r1-d,ln);	% Coordonn�es L1 et L2 dans le rep�re normalis�
	L2=dot(r2-d,ln);

    % Calcul de G0 et H0 (c'est compact, il faut r�f�rer � mes notes pour les d�tails)
	%
	% Formules:
	%
	% G0=1/(2*pi)*( atan(L2/D) - atan(L1/D) )
	% H0=1/(2*pi)*( L2*(log(|r2|)-1) - L1*(log(|r1|)-1) + D*(atan(L2/D) - atan(L1/D) )
	%
	G0=0;
    H0=0;
    if abs(D)>seuil
		T2=atan(L2/D);
        T1=atan(L1/D);
        H0=(T2-T1);
		G0=D*H0;
    else
        D=0;
    end
	if abs(L2)>seuil
		G0=G0+L2*(log(norm(r2))-1);
	end
	if abs(L1)>seuil
		G0=G0-L1*(log(norm(r1))-1);
	end
	G0v(i)=G0; 
    H0v(i)=H0;
end

% Facteur 2*pi pour tous les coefficients calcul�s (1 par �l�ment)
G0v=G0v/(2*pi); 
H0v=H0v/(2*pi);
