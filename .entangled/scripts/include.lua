function Para(elem)
    if #elem.content < 3 then
        return
    end
    if elem.content[1].t ~= "Str" or elem.content[1].text ~= "!include" then
        return
    end
    if elem.content[3].t ~= "Str" then
        print("Don't know how to include: "..elem.content[3])
        return
    end
    filename = elem.content[3].text
    content = io.input(filename):read("a")
    if filename:sub(-2) == "md" then
        return pandoc.read(content).blocks
    end
    if filename:sub(-4) == "html" then
        return pandoc.RawBlock("html", content)
    end
end
