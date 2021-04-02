
function popup_window(modal_ID, action)
{
	/* tool function for pop up views */
	// modal_ID: pop up content ID
	// action: open or close
	
	var	modal = document.getElementById(modal_ID);
	
	if(action == "open"){
		modal.style.display = "block";
		
	}else if(action == "close"){
		modal.style.display = "none";
		
	}	
}


function alert_checkbox_select(checkbox_ID, message)
{
	/* reverse the decision if the user is not sure */
	// checkbox_ID: checkbox clicked
	// message: warning message to display
	
	var	checkbox = document.getElementById(checkbox_ID);
	
	if(checkbox.checked && !confirm(message))
	{
		checkbox.checked = false;
	}
}


function click_show_hide(checkbox_ID, target_ID)
{
	/* show or hide a block controlled by a checkbox */
	// checkbox_ID: control checkbox
	// target_ID: target block
	
	var	checkbox = document.getElementById(checkbox_ID),
		target = document.getElementById(target_ID);
	
	if(checkbox.checked){
		target.style.display="block";
	}else{
		target.style.display="none";
	}
}


// Two functions for highlighting content with mouse over and out
function highlight_mouseover(target)
{
	var i, targets = target.split(',')
	
	for(i=0; i<targets.length; i++)
	{
		document.getElementById(targets[i]).style.color = "red";
	}
}

function highlight_mouseout(target)
{
	var i, targets = target.split(',')
	
	for(i=0; i<targets.length; i++)
	{
		document.getElementById(targets[i]).style.color = "black";
	}
}
