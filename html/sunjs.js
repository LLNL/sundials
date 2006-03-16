/*
 *******************************************************
 Script to expand/contract announcements table
 *******************************************************
 */ 

function doTable(i){
  if(!document.getElementById)
      return;

  var d=document.getElementById(i);

  if(d.style.display=='none')
  {
   d.style.display='block';
   document ['arrowImg'].src = "img/arrowUp.png";
  }
  else
  {
   d.style.display='none';
   document ['arrowImg'].src = "img/arrowDown.png";
  }
}

/*
 *******************************************************
 Script to change images on rollover
 *******************************************************
 */ 

function img_act(imgName) {
      imgOn = eval(imgName + "on.src");
      document [imgName].src = imgOn;
}
function img_inact(imgName) {
      imgOff = eval(imgName + "off.src");
      document [imgName].src = imgOff;
}

/*
 *******************************************************
 Script to read a cookie and redirect based on its value 
 *******************************************************
 */ 

function createCookie(name,value,days)
{
	if (days)
	{
		var date = new Date();
		date.setTime(date.getTime()+(days*24*60*60*1000));
		var expires = "; expires="+date.toGMTString();
	}
	else var expires = "";
	document.cookie = name+"="+value+expires+"; path=/";
}

function readCookie(name)
{
	var nameEQ = name + "=";
	var ca = document.cookie.split(';');
	for(var i=0;i < ca.length;i++)
	{
		var c = ca[i];
		while (c.charAt(0)==' ') c = c.substring(1,c.length);
		if (c.indexOf(nameEQ) == 0) return c.substring(nameEQ.length,c.length);
	}
	return null;
}

function eraseCookie(name)
{
	createCookie(name,"",-1);
}

function testRegistration(which, ver)
{
        if (readCookie('SunRegCookie') == 'OK')
        {
            if (ver == 'NEW')
            { 
                var dwnld_file = "code/" + which +".tar.gz"; 
                location.href = dwnld_file;
            }
            else
            {
                var dwnld_file = "code/old/" + which +".tar.gz"; 
                location.href = dwnld_file;
            }
        }
        else
        {
                createCookie('SunPckgCookie',which);
                createCookie('SunVersionCookie',ver);
                location.href = 'sundials_agree.html';
        }

}

function setFormValues()
{
        var package_name = readCookie('SunPckgCookie');
        var package_ver  = readCookie('SunVersionCookie');
        if (package_ver == 'NEW')
        {
                document.sunRegForm.packageName.value = package_name;
                document.sunRegForm.subject.value = "SUNDIALS Download Registration (" + package_name + ")";
                document.sunRegForm.redirect.value = "../download/download.html";
        }
        else
        {
                document.sunRegForm.packageName.value = package_name;
                document.sunRegForm.subject.value = "SUNDIALS Download Registration - OLD version (" + package_name + ")";
                document.sunRegForm.redirect.value = "../download/oldversions.html";
        }

}