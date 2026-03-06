function dydt = oregonator_rhs(~,y,p)
%% Paper Title: "An Efficient Adjoint-Free Melnikov-Lyapunov Diagnostic for
%%               Transition Detection in Slow-Fast Chemical Reaction Networks"
%% Author 1:    Sri Venkata Durga Sudarsan Madhyannapu
%% Author 2:    Pradheep Kumar S.
%%
%% Affiliation 1: Freshmen Engineering Department, NRI Institute of Technology
%%                (Autonomous), Pothavarappadu, Agiripalli, Vijayawada,
%%                521212, Andhra Pradesh, India
%% Affiliation 2: Research Scholar, Jawaharlal Nehru Technological University
%%                Kakinada, Andhra Pradesh, India
%% Affiliation 3: School of Basic Sciences, SRM University AP, Neerukonda,
%%                Mangalagiri Mandal, Guntur, 522240, Andhra Pradesh, India
%%
%% Journal:       Communications in Nonlinear Science and Numerical Simulation
%%                Elsevier, ISSN: 1007-5704
%% Manuscript ID: CNSNS-D-26-00848
%% Submitted:      19 February 2026
%% Status:        Under Review, 2026
%% SSRN:          https://ssrn.com/abstract=6275667
%% SSRN ID:        6275667 (Distributed: 02/20/2026)
%%
%% Paper Model 6: Modified Oregonator (BZ reaction)
%% Parameters: q=0.0008, f=1.4, w=0.1, eps=4e-4, delta=2e-4
if nargin<3, p.eps=4e-4;p.delta=2e-4;p.q=0.0008;p.f=1.4;p.w=0.1;p.mu=0; end
dydt=zeros(3,1);
dydt(1)=(1/p.eps)*(y(2)-y(1)*y(2)+y(1)-p.q*y(1)^2)+p.mu*y(1);
dydt(2)=(1/p.delta)*(-y(2)-y(1)*y(2)+p.f*y(3));
dydt(3)=p.w*(y(1)-y(3));
end
