Title = "unknown"

function Meta(meta)
    Title = meta.subtitle
end

function Pandoc(doc)
    return pandoc.Pandoc{Title}
end
