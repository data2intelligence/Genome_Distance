
//Grab and format a column from the table
function get_table_cell(j, col)
{
	var $col = $(col), $text = '', input = $col.find('input'), select = $col.find('select');
	
	
	if(input.length > 0)
	{
		// divide by input field types	
		if (input[0].type == 'text'){
			$text = input[0].value.trim().replace('"', '""');
			
		}else if(input[0].type == 'radio' || input[0].type == 'checkbox'){
			$text = input[0].checked;
			
		}else{
			$text = null;
		}
		
		
	}else if (select.length > 0){
		// option selection 
		$text = select[0].options[select[0].selectedIndex].value;
		
	}else{
		// regular fixed cell content
		$text = $col.text().trim().replace('"', '""');
	}
	
	return $text; // escape double quotes
}


function export_table_csv($table, filename)
{
	var $headers = $table.find('tr:has(th)'), $rows = $table.find('tr:has(td)'),
	// Temporary delimiter characters unlikely to be typed by keyboard
	// This is to avoid accidentally splitting the actual contents
	tmpColDelim = String.fromCharCode(11), // vertical tab character
	tmpRowDelim = String.fromCharCode(0), // null character

	// actual delimiter characters for CSV format
	colDelim = '","', rowDelim = '"\n"';

	// Grab text from table into CSV formatted string
	var csv = '"';
	
	csv += formatRows($headers.map(grabRow));
	csv += rowDelim;
	
	var row_content = $rows.map(grabRow);
	
	// remove empty rows
	row_content = row_content.map(
			function(i, v){return v.trim() == "" ? null:v;}
			);
	
	csv += formatRows(row_content) + '"';

	// Data URI
	var csvData = 'data:application/csv;charset=utf-8,' + encodeURIComponent(csv);

	// For IE (tested 10+)
	if (window.navigator.msSaveOrOpenBlob) {
		var blob = new Blob([ decodeURIComponent(encodeURI(csv)) ], {
			type : "text/csv;charset=utf-8;"
		});
		navigator.msSaveBlob(blob, filename);
	} else {
		// ,'target' : '_blank' //if you want it to open in a new window
		$(this).attr({
			'download' : filename,
			'href' : csvData
		});
	}
	
	//------------------------------------------------------------
	// Helper Functions
	// ------------------------------------------------------------
	// Format the output so it has the appropriate delimiters
	function formatRows(rows) {
		return rows.get().join(tmpRowDelim).split(tmpRowDelim).join(rowDelim)
				.split(tmpColDelim).join(colDelim);
	}

	// Grab and format a row from the table
	function grabRow(i, row) {
		var $row = $(row);

		// for some reason $cols = $row.find('td') || $row.find('th') won't work...
		var $cols = $row.find('td');

		if (!$cols.length)
			$cols = $row.find('th');

		return $cols.map(get_table_cell).get().join(tmpColDelim);
	}
}



function table_search_filter(document, ID_search_filter_keyword, ID_select_filter_column, ID_table, flag_regex)
{
	/* filter table rows by search */
	
	var	table = document.getElementById(ID_table),
		filter_pattern = document.getElementById(ID_search_filter_keyword).value.trim(),
		target = document.getElementById(ID_select_filter_column),
		tr = table.getElementsByTagName("tr"),
		td, i, txtValue;
	
	if (!flag_regex){
		filter_pattern = escapeRegExp(filter_pattern);
	}
	
	// convert to regular expression
	filter_pattern = new RegExp(filter_pattern, "i");
	
	if(target.selectedIndex == 0){
		window.alert("Please select a target column for the search filter!")
	
	}else{
		target = target.options[target.selectedIndex].value;
		
		// Loop through all table rows, and hide those who don't match the search query
		for (i = 0; i < tr.length; i++)
		{
			td = tr[i].getElementsByTagName("td")[target];
			
			if (td){
				//txtValue = td.textContent || td.innerText;
				txtValue = get_table_cell(target, td);
				
				// 
				// txtValue.toUpperCase().indexOf(filter) > -1
				if (txtValue.match(filter_pattern) != null) {
					tr[i].style.display = "";
				} else {
					tr[i].style.display = "none";
				}
			}
		}
	}
}



function delete_table_column(table_ID, select_delete_column)
{
	/* A tool function to delete table column within the HTML table */
	// table: table ID
	// select_delete_column: option selection ID for target column 
	
	var select_col = document.getElementById(select_delete_column);
	
	if(select_col == null){
		window.alert("Error: No selection ID " + select_delete_column);
		
	}else if(select_col.selectedIndex <= 0){
		window.alert("Error: No column to delete!");
	
	}else{
		// Remove tablesorter from the table before updating the header, trigger destroy will restore the header setting to an initial state
		$('#' + table_ID).trigger('destroy'); 
		
		select_col = select_col.options[select_col.selectedIndex].value;
		
		$('#' + table_ID + ' tr').find('td:eq(' + select_col + '),th:eq(' + select_col + ')').remove();
		
		load_table_columns_to_options();
	}
}

function append_table_column(table_ID, column_append_name)
{
	/* A tool function to append table column within the HTML table */
	// table: table ID
	// column_append_name: option selection ID for column name 
	
	var append_col = document.getElementById(column_append_name);
	
	if(append_col == null){
		window.alert("Error: No selection ID " + column_append_name);
		
	}else{
		append_col = append_col.value.trim();
		
		if(append_col.length == 0){
			window.alert("Please input some column name!");
			
		}else{				
			// Remove tablesorter from the table before updating the header, trigger destroy will restore the header setting to an initial state 
			$('#' + table_ID).trigger('destroy'); 
			
			$('#' + table_ID + ' tr').each(function(i){
				if(i==0){
					$(this).append("<th class='input'>" + append_col + "</th>");
				}else{
					$(this).append("<td><input type='text'></td>");
				}
			});
			
			load_table_columns_to_options();
		}
	}
}

function manipulate_table_columns(table, select_src, select_dst, input_separator, action, text = null, target = null, vocabulary_map=null)
{
	/* A tool function to manipulate table columns within the HTML table */
	// table: table ID
	// select_src: option selection ID for source column
	// select_dst: option selection ID for destination column
	// input_separator: input text ID for separator when joining column text together
	// action: type of actions: copy, concatenate, append, clear, separate, text (input text directly), replace (replace by regex)
	// text input: if any, such as from clipboard, or regex for replace
	// target: if any, regex replace target
	
	var	src = $('#' + select_src)[0],
		dst = $('#' + select_dst)[0],
		separator = $('#' + input_separator)[0].value,
		start_state = $('#' + table).html(), i, fields;
	
	if(src.selectedIndex == 0 && ["copy", "concatenate"].includes(action)){
		window.alert("Please select a source column!")
		
	}else if(dst.selectedIndex == 0){
		window.alert("Please select a destination column!")
		
	}else{
		src = src.options[src.selectedIndex].value;
		dst = dst.options[dst.selectedIndex].value;
		
		if(text != null && action == "paste"){
			text = text.split(String.fromCharCode(13));
		}	
		
		// Copy column cells
		$('#' + table + ' tbody tr').each(function(i){
			
			// jump filtered rows
			if($(this).attr('class') == 'filtered') return true;
			
			var input_box = $(this).children()[dst].getElementsByTagName('input')[0];
			
			if(action == "copy"){
				input_box.value = $(this).children()[src].innerText;
				
			}else if(action == "concatenate_front"){
				input_box.value = $(this).children()[src].innerText + separator + input_box.value;
				
			}else if(action == "concatenate_rear"){
				input_box.value += separator + $(this).children()[src].innerText;
				
			}else if(action == "append_front"){
				input_box.value = separator + input_box.value;
				
			}else if(action == "append_rear"){
				input_box.value += separator;
				
			}else if(action == "clear"){
				input_box.value = '';
				
			}else if(action == "paste"){
				if(text != null)
				{
					if(text.length > 1 && i < text.length){
						input_box.value = text[i].trim();
						
					}else if (text.length == 1){
						input_box.value = text[0].trim();
					}
				}
			
			}else if(action == "replace"){
				input_box.value = input_box.value.replace(text, target).trim();
				
			}else if(action == "separate_rear"){
				fields = input_box.value.split(text);
				input_box.value = fields[fields.length - 1].trim();
				
			}else if(action == "separate_front"){
				input_box.value = input_box.value.split(text)[0].trim();
				
			}else if(action == "translate" && vocabulary_map !=null){
				for(i=0; i < vocabulary_map.length; i++)
				{
					if(input_box.value.match(vocabulary_map[i][0]) != null)
					{
						input_box.value = input_box.value.replace(vocabulary_map[i][0], vocabulary_map[i][1]).trim();
					}
				}
				
			}
			
		});	
		
		// this is needed for table sorter
		$("#" + table).trigger("update");
	}
}


// Table Undo Redo Stack
function UndoStack_Table(table)
{
	this.table = table;	// ID of current table
	this.stack = [];	// clones of each states
	this.current = -1;	// index to current element
};


UndoStack_Table.prototype.push = function ()
{
	this.current++;	// wind to next position
	
	this.stack.splice(this.current);
	
	// update to enable a normal html() call later
	$('#' + this.table).find("input").each(
		function(){
			if (this.defaultValue !== this.value){this.defaultValue = this.value;}
		});
	
	this.stack.push($('#' + this.table).html());
};


UndoStack_Table.prototype.undo = function ()
{
	// the first position is already the initial elements
	if (this.current > 0){
		this.current--;
		var item = this.stack[this.current];
		
		// have to destroy tablesorter to regenerate new one later
		$('#' + this.table).trigger('destroy');
		$('#' + this.table).html(item);
		load_table_columns_to_options();
	}else{
		alert("Already at oldest change");
	}
};


UndoStack_Table.prototype.redo = function ()
{
	this.current++;
	var item = this.stack[this.current];
	
	if (item) {
		// have to destroy tablesorter to regenerate new one later
		$('#' + this.table).trigger('destroy');
		$('#' + this.table).html(item);
		load_table_columns_to_options();
	}else{
		alert("Already at newest change");
	}
};
