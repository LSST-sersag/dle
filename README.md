####  2021 Enabling Science Call project 
# Building Deep Learning Engine (DLE) for AGN light-curves

![image001](https://user-images.githubusercontent.com/78701856/121324461-0b87ac80-c911-11eb-8196-2c688f61a4bd.png)

<b> PIs:</b>  Andjelka Kovacevic (andjelka@math.rs), Dragana Ilic (dilic@math.rs), University of Belgrade,

<b> Co-Is:</b>  Luka C. Popovic (Astronomical Observatory Belgrade, University of Belgrade, lpopovic@aob.rs), Paula Sánchez Sáez (Pontificia Universidad Católica de Chile, pasanchezsaez@gmail.com), Robert Nikkuta
(NOIRLab, robert.nikutta@noirlab.edu)

<b> Postdocs: </b> Viktor Radović (rviktor@math.rs,, University of Belgrade), Djordje Savić(djasvic@aob.rs,,
Astronomical Observatory Belgrade)

<b> PhD students: </b> Iva Čvorović Hajdinjak(iva_cvorovic@gmail.com,University of Belgrade), Isidora Jankov (isidora_jankov@math.rs,University of Belgrade), Nemanja Rakić (rakinemanja@gmail.com, University of Belgrade)

<b> Summary: </b> <i>  Developing deep learning engines (DLEs) for non-parametric modeling and extracting of information from active galactic nuclei (AGN) light-curves (LCs), which are directly related to the scientific objectives of the LSST Exploring transient optical sky. Developed DLEs Jupyter notebooks might be adaptable for modeling of light-curves of other objects, and will be tested on LSST Data Previews, as well as on other datasets through NOIRLab services. </i>


<b>  PROJECT DESCRIPTION, MOTIVATION, AND RELEVANCE TO LSST  </b>

<b> Motivation and relevance. </b> Exploring transient optical sky (ETOS) is among four Rubin Observatory Legacy Survey of Space and Time (LSST) key science drivers (Ivezić et al. 2019). We are motivated by the ETOS LSST science opportunity No. 14 (see item 14 in section 4, Ivezić et al. 2019), which enables to harness LSST light curves (LC) of active galactic nuclei (AGN) for photometric reverberation mapping (PhotoRM, Chelouche & Daniel 2012). Our goal is to build a deep learning engine (DLE) for LC nonparametric modeling and implementation of the PhotoRM procedure to respond to the observing strategy of the LSST (Jones et al. 2020). Our student has already started experiments of DLE LC modeling, and pilot results are shown in Figure 1. Here the DLE has been trained on object 1H0419m577 LC and produces the predicted points in green. Learned LC will
enable us to improve time-lag determination as a goal of
PhotoRM, but also to hunt for Transient Fueling Events:
Temporary AGNs and Cataclysmic AGN Outbursts (see pp.
369-371 in LSST Collaboration Science book 2009). Further,
our group has already begun PhotoRM calculations and
some initial procedures have been tested in Kovacevic et al.
(submitted 2021). Our DLE will be developed as user-
friendly Jupyter notebooks, in such a way to be adaptable to
other object LC, and tested on LSST Data Previews (e.g.,
DP0). Initial stages of our DLE will be presented as an
ePoster  at the  European Astronomical Society  Annual
Meeting - EAS 2021, subtheme SS32: Machine learning and visualisation in data intensive era (Jun 28 - Jul 2, 2021), and at the 13th SCSLSA Conference (Aug 23-27, 2021). We plan to further disseminate results at international conferences and in peer-reviewed publications.

<img width="473" alt="image" src="https://user-images.githubusercontent.com/78701856/121322214-14777e80-c90f-11eb-95d1-bd5e2e3a0c80.png">

<b> Project description. </b> For the DLEs developing and testing, we propose the following: (i) 2 long-term projects for already enrolled PhD students lasting for the duration of the grant period (~10 months), (ii) 4 short-term internships for undergraduate/master students (~3 months each) starting mid- project, and (iii) a series of tutorials/hack days/workshops, with special focus on local and regional involvement.
2x Long-term student projects description:
<i> PhD student DLE1. </i> Non-parametric modelling of AGN LC across redshifts: using different techniques: Gaussian processes, artificial neural networks. DLE environment in Jupyter notebooks through NOIRLab services. Contribution of notebooks to the Astro Data Lab science platform.


<i> PhD student DLE2. </i> PhotoRM in action. Development of Jupyter notebooks for extraction of time-lags from photometric AGN light curves across redshift based on pure PRM principles. Contribution of notebooks to the Astro Data Lab science platform.
In the first phase, DLEs will be developed on the local Cluster. We will run and test computationally intensive codes on cluster GPUs, for which we plan an add-on through this project. DLEs will be tested on available datasets (ZTF, CRTS) through NOIRLab services, and on the LSST Data Previews and artificially generated LC (such as PLAsTiCC 2). Since our LC modeling is data driven (non-parametric) both projects can be tested on other objects and transients LC, such as binary stars and microlensing. Trimester progress reports and the summary reports will be provided.

<i> 4 x Short-term student internships. </i> Four undergraduate student internships for user testing of DLE1 and DLE2 projects on different datasets (ZTF, CRTS) and show-cases, providing solutions for minor problems and bugs. New students from physics, astrophysics, informatics, and statistics background will be selected based on their motivation and skills from the PIs institutions, and the region. Internship will start mid-project, between university semesters. The expected internship duration is 13 weeks each, yielding approx. 200 working hours per student, and will be funded with $US 500/internship. Funded students are expected to submit a summary report within a month of completing the project.

<b> Supervision. </b> Two students (1 PhD, 1 undergraduate) will be supervised by each of the two PIs, with the support of other co-PIs. The supervisor’s role is to formulate their research project and related questions, decide what methods of research to use, evaluate the results of their research, ensure their work meets the necessary standards expected by Rubin LSST Community and academia, help to meet the deadlines, use feedback to enhance their work, and provide guidance for preparation of conference presentations and journal papers. The supervision will also teach shared responsibilities and mutual respect. 

<b> HR and Management.</b>  Our team consists of five experts in AGN time-domain studies supporting diversity in any sense, gender, racial, national, etc., emphasizing female participation. The team has successfully managed many scientific projects and has a long and rich experience in organizing events of various kinds and formats.


<b> Deliverables: </b>
<ol> <li> <b> Presentations at conferences: </b>
1xDLE1, at the  European Astronomical Society  Annual Meeting - EAS 2021, subtheme SS32: Machine
learning and visualisation in data intensive era (Jun 28 - Jul 2, 2021) 1xDLE1, 1XDLE2 at the 13th SCSLSA Conference (Aug 23-27, 2021)
1xDLE1, 1xDLE2 in 2022 (SPIG 2022, LSST PCW 2022;
  </li>
  <li>
    <b> Papers: </b>1xDLE1, 1xDLE2 in Astronomische Nachrichten by the end of 2021; 1xDLE1, 1xDLE2 in a scientific journal by the end of 2022; </li>
  <li> Open-source Jupyter notebooks: DLE1, DLE2 by May 2022. </li>


  <li>  <b> Tutorials/Hack-Days. </b>  A series of online tutorials/hack-days will be organized regularly. We envision at least three (but preferentially more) online events every trimester. They will consist of two interconnected modules, one to work together on developing and solving required tasks, and the other on guide-through tutorials and showcases. The focus is on engaging the local and regional participation, but international participants are welcome. The tutorials will use supporting online tools and platforms, with the focus on NOIRLab data services. The key presenters will be students, and the sessions will provide them space to develop and practice their communication and management skills. Hack-days reviewers will be external software experts. </li>
  
  <li> <b> Final workshop. </b> Towards the end of the project (~ in 10th month) a final 2-day online workshop will be organized with the aim to review, evaluate, and disseminate the developed tool(s). Event will have a format of 1-2 invited talks and short student presentation, with a wide participation from the Rubin LSST Community. </li>
  </ol>


<b> Team members </b> 

<img align = "center" width="840" alt="GroupPhoto" src="https://user-images.githubusercontent.com/78701856/121871316-41a7a080-cd04-11eb-9bcf-d34f6482cebc.png">

<b> Sponsors </b>
<br>
This project is graciously supported by a grant from the LSST Corporation's Enabling Science Small Grants Program and by the University of Belgrade.



<table>
  
  <tr>
    <td align="center"><img width=220 height=100 src="https://user-images.githubusercontent.com/78701856/121324461-0b87ac80-c911-11eb-8196-2c688f61a4bd.png"> </td>
    <td align="center"> <img width =210 height=100 src="https://user-images.githubusercontent.com/78701856/121673183-edf14900-cab0-11eb-83e6-ceeac881c2f9.png"></td> 
    <td align="center"> <img width=220 height=100 src="https://user-images.githubusercontent.com/78701856/121673230-f8abde00-cab0-11eb-8f1c-eacf6c399c4b.png"> </td>
  </tr>
  <tr>
    <td align="center">  <img width=250 height=100 src="https://user-images.githubusercontent.com/78701856/121677377-116ac280-cab6-11eb-80a9-058141fd9939.png"> </td>
    <td align="center"> <img width=110 height=100 src="https://user-images.githubusercontent.com/78701856/121677418-1d568480-cab6-11eb-8c42-0a72b4e8ead5.gif"> </td>
    <td align="center"> <img width=200 height=100 src="https://user-images.githubusercontent.com/78701856/121677472-3101eb00-cab6-11eb-87fc-b1c193d03ef0.png"> </td>

  </tr>
 </table>


