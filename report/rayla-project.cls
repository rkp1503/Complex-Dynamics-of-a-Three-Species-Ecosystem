%-----------------------------------------------------------------------------
% Author: Rayla Kurosaki
% GitHub: https://github.com/rkp1503
%-----------------------------------------------------------------------------
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{rayla-project}[Rayla's Project's Class]
\LoadClass{article}

%-----------------------------------------------------------------------------
% 3rd party packages
%-----------------------------------------------------------------------------
\usepackage{amsmath}
\usepackage[backend=biber,style=numeric,sorting=ynt]{biblatex}
\usepackage{geometry}
\usepackage{graphicx,epstopdf}
\usepackage{hyperref}
\usepackage{ifthen}
\usepackage{lipsum}

%-----------------------------------------------------------------------------
% Geometry
%-----------------------------------------------------------------------------
\geometry{left=2cm,right=2cm,top=2cm,bottom=2cm}
\setlength\parindent{0pt}

%-----------------------------------------------------------------------------
% Setting figures directory
%-----------------------------------------------------------------------------
\graphicspath{{figures/}}

%-----------------------------------------------------------------------------
% Custom Title Style
%-----------------------------------------------------------------------------
\newcommand{\myUniversity}[1]{\newcommand{\myUniversityLoc}{#1}}
\newcommand{\myTitle}[1]{\newcommand{\myTitleLoc}{#1}}
\newcommand{\myName}[1]{\newcommand{\myNameLoc}{#1}}
\newcommand{\courseid}[1]{\newcommand{\courseidLoc}{#1}}
\newcommand{\courseName}[1]{\newcommand{\courseNameLoc}{#1}}
\newcommand{\professorName}[1]{\newcommand{\professorNameLoc}{#1}}
\newcommand{\term}[1]{\newcommand{\termLoc}{#1}}
\newcommand{\dueDate}[1]{\newcommand{\dueDateLoc}{#1}}
    
\renewcommand{\maketitle}{
    \thispagestyle{plain}
    \begin{center}
        {\scshape \myUniversityLoc} \\[10pt]
        \hrule
        \vspace{10pt}
        {\LARGE \bfseries \myTitleLoc} \\[5pt]
        Student: {\itshape \myNameLoc}
        \vspace{8pt}
        \hrule
        \vspace{10pt}
        {Course: {\itshape \courseidLoc \ (\courseNameLoc)}} \hfill {Term: {\itshape \termLoc}}\\
        {Professor: {\itshape\professorNameLoc}} \hfill {Due date: {\itshape\dueDateLoc}}
        \vspace{20pt}
    \end{center}
}

%-----------------------------------------------------------------------------
% Custom referencing
%-----------------------------------------------------------------------------
\DeclareRobustCommand{\myref}[2][]{%
    \ifthenelse{\equal{#1}{Equation} \OR \equal{#1}{System} \OR \equal{#1}{Model} \OR \equal{#1}{Matrix}}
    {#1~\eqref{#2}}
    {#1~\ref{#2}}%
}

%-----------------------------------------------------------------------------
% Inclusion of section number in equation numbering
%-----------------------------------------------------------------------------
\numberwithin{equation}{section}

%-----------------------------------------------------------------------------
% Bibliography settings
%-----------------------------------------------------------------------------
\addbibresource{main.bib}
