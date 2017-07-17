"""
Description:
    what does this script module do? How to do it.
    Hi Fei

    This is how I calculate bostick resistivity and depth
Author: fei.zhang@ga.gov.au

Date: 2017-07-14
"""



	 private double bostick_depth(double f, double rho){
		 double mu = 4*Math.PI*Math.pow(10,-7);
		 double h = Math.sqrt(rho/(f*2*Math.PI*mu));
		 //System.out.println(h);
		 return h;
	 }


private double bostick__resistivity(double f, double rho,double pha_radian){
		 // phase in radian
		 //if (pha_radian > 180)
		 pha_radian = pha_radian*Math.PI/180.0;
		 double rho_b = rho*(Math.PI/(2*pha_radian) -1);
		 //System.out.println(rho_b);
		 return rho_b;
	 }



-----Original Appointment-----
From: Zhang Fei
Sent: Thursday, 13 July 2017 4:25 PM
To: Wang Liejun
Subject: chat about DOI [SEC=UNCLASSIFIED]
When: Friday, 14 July 2017 11:00 AM-12:00 PM (UTC+10:00) Canberra, Melbourne, Sydney.
Where: Fei's work area



DOI and
ModEM input
other MT things.
