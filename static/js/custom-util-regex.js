// Regular expression tools

function get_regex_lst(keywords)
{
	var regex_lst = keywords.split(',');
	
	for(var i=0; i< regex_lst.length; i++)
	{
		regex_lst[i] = new RegExp(regex_lst[i].trim(), "gi"); // global match, ignore case	
	}
	
	return regex_lst;
}


function highlight_words(regex_lst, element)
{
	/* high light text */
	// regex_lst: pattern list to search
	// element: target to check
	
	$(regex_lst).each(function(i) {
		regex = this;
		
		var textNodes = $(element).contents().filter(function() { return this.nodeType === 3 });
				
		textNodes.each(function(j) {
			var content = $(this).text();
			content = content.replace(regex, '<span class="highlight">' + '$&' + '</span>');
			$(this).replaceWith(content);
		});
	});
	
}

function escapeRegExp(text)
{
	// Escape all special characters in regular expression for string absolute match mode
	return text.replace(/[-[\]{}()*+?.,\\^$|#\s]/g, '\\$&');
}


function convertRegExp(pattern)
{
	// convert strings with /text/flag to a regular expression, or still string exact if not
	if(pattern[0] == '/'){
		pattern = pattern.split('/');
		
		if(pattern.length >= 3)
		{
			flag = pattern[pattern.length-1];
			pattern = pattern.slice(1,-1).join('/');
			
			return new RegExp(pattern, flag);
		}
	}
	
	// exact match
	pattern = escapeRegExp(pattern)
	
	return new RegExp(pattern, 'g');
}
