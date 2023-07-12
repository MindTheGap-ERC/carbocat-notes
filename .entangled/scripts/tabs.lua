function Div(elem)
    if elem.classes[1] ~= "tabs" then
        return
    end
    local buttons = {}
    local tabs = {
    }

    elem:walk({
        Div = function (elem)
            if elem.classes[1] ~= "tab" then
                return
            end
            if elem.identifier == nil then
                print("This tab has no id.")
                print(elem)
                return
            end
            local title = elem.attr["title"]
            if title == nil then
                title = elem.identifier
            end
            table.insert(buttons, pandoc.RawBlock("html",
                "<button id=\"".. elem.identifier .."-button\" class=\"tablink\" onclick=\"openTab(event, '".. elem.identifier  .."')\">".. 
                title .."</button>"))
            table.insert(tabs, elem)
        end
    })

    table.insert(tabs, pandoc.RawBlock("html", "<script>"..
        "document.getElementById(\"".. tabs[1].identifier .."-button\").click();"..
        "</script>"))
    table.insert(tabs, 1, pandoc.Div(buttons, {class="buttons"}))
    return pandoc.Div(tabs, {class="tabs"})
end

