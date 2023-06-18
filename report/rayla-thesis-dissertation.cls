%-----------------------------------------------------------------------------
% Author: Ramsey (Rayla) Phuc
% Alias: Rayla Kurosaki
% GitHub: https://github.com/rkp1503
%-----------------------------------------------------------------------------
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{rayla-thesis-dissertation}[Rayla's Thesis and Dissertation Class]
\LoadClass[twoside]{report}

%-----------------------------------------------------------------------------
% 3rd party packages
%-----------------------------------------------------------------------------
\usepackage{amsmath, amsthm, mathtools}
\usepackage[backend=biber,style=numeric,sorting=nyt]{biblatex}
\usepackage{fancyhdr}
\usepackage{float}
\usepackage{geometry}
\usepackage{graphicx,epstopdf}
\usepackage{hyperref}
\usepackage{ifthen}
\usepackage{lipsum}
\usepackage{parskip}

%-----------------------------------------------------------------------------
% Defining new commands
%-----------------------------------------------------------------------------
\newcommand{\myUniversity}[1]{\newcommand{\myUniversityLoc}{#1}}
\newcommand{\myCollege}[1]{\newcommand{\myCollegeLoc}{#1}}
\newcommand{\myDepartment}[1]{\newcommand{\myDepartmentLoc}{#1}}
\newcommand{\myLogo}[1]{\newcommand{\myLogoLoc}{#1}}
\newcommand{\myTitle}[1]{\newcommand{\myTitleLoc}{#1}}
\newcommand{\myName}[1]{\newcommand{\myNameLoc}{#1}}
\newcommand{\myAdvisor}[1]{\newcommand{\myAdvisorLoc}{#1}}
\newcommand{\submissionDate}[1]{\newcommand{\submissionDateLoc}{#1}}

%-----------------------------------------------------------------------------
% Geometry
%-----------------------------------------------------------------------------
\geometry{left=2cm,right=2cm,top=2cm,bottom=2cm}
\setlength\parindent{0pt}
\allowdisplaybreaks

%-----------------------------------------------------------------------------
% Setting figures directory
%-----------------------------------------------------------------------------
\graphicspath{{figures/}}
    
%-----------------------------------------------------------------------------
% Custom Title Style
%-----------------------------------------------------------------------------
\renewcommand{\maketitle}{
    \thispagestyle{empty}
    \begin{center}
        \newcommand{\HRule}{\rule{\linewidth}{0.5mm}}
        \textsc{\LARGE \myUniversityLoc}\\[0.5cm]
        \textsc{\Large \myCollegeLoc}\\[0.5cm]
        \textsc{\large \myDepartmentLoc}\\[0.5cm]
        \begin{figure}[H]
            \centering
            \includegraphics[scale=0.15]{\myLogoLoc}
            \label{fig:logo}
        \end{figure}
        \HRule\\[0.4cm]
        {\huge \textbf{\myTitleLoc}}\\[0.4cm]
        \HRule\\[1.5cm]
    \end{center}
    \begin{align*}
        &\text{By:} & &\text{\myNameLoc}\\
        &\text{Research Advisor:} & &\text{\myAdvisorLoc}\\
        &\text{Submission Date:} & &\text{\submissionDateLoc}
    \end{align*}
}

%-----------------------------------------------------------------------------
% Custom Headers and Footers
%-----------------------------------------------------------------------------
\pagestyle{fancy}
\renewcommand{\headrulewidth}{0.4pt}
\fancyhead{}
\fancyhead[LE]{\myTitleLoc}
\fancyhead[CE]{}
\fancyhead[RE]{}
\fancyhead[LO]{}
\fancyhead[CO]{}
\fancyhead[RO]{\myTitleLoc}
\renewcommand{\footrulewidth}{0.4pt}
\fancyfoot{}
\fancyfoot[LE]{\thepage}
\fancyfoot[CE]{}
\fancyfoot[RE]{\myNameLoc}
\fancyfoot[LO]{\myNameLoc}
\fancyfoot[CO]{}
\fancyfoot[RO]{\thepage}

%-----------------------------------------------------------------------------
% Custom referencing
%-----------------------------------------------------------------------------
\DeclareRobustCommand{\myref}[2][]{%
    \ifthenelse{\equal{#1}{Equation} \OR \equal{#1}{System} \OR \equal{#1}{Model} \OR \equal{#1}{Matrix} \OR \equal{#1}{Theorem}}
    {#1~\eqref{#2}}
    {#1~\ref{#2}}%
}


% Theorem Styles
\theoremstyle{plain}
\newtheorem{theorem}{Theorem}
\newtheorem{corollary}{Corollary}[theorem]
\newtheorem*{main}{Main Theorem}
\newtheorem{lemma}{Lemma}
\newtheorem{claim}{Claim}
\newtheorem{proposition}{Proposition}
\theoremstyle{definition}
\newtheorem{define}{Definition}
\theoremstyle{remark}
\newtheorem*{notation}{Notation}
\theoremstyle{remark}
\newtheorem*{remark}{Remark}
\theoremstyle{remark}
\newtheorem{example}{\textbf{Example}}%[section]

%-----------------------------------------------------------------------------
% Inclusion of chapter number in equation numbering
%-----------------------------------------------------------------------------
\numberwithin{equation}{chapter}
\numberwithin{theorem}{chapter}

%-----------------------------------------------------------------------------
% Bibliography settings
%-----------------------------------------------------------------------------
\addbibresource{main.bib}
