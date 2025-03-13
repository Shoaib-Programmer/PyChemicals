// @ts-check 
import { defineConfig } from "astro/config"; 
import starlight from "@astrojs/starlight"; 
 
export default defineConfig({ 
  image: { 
    service: { entrypoint: "astro/assets/services/noop" }, 
  }, 
  site: "https://Shoaib-Programmer.github.io/PyChemicals", 
  base: "/PyChemicals/", 
  integrations: [ 
    starlight({ 
      title: "PyChemicals", 
      logo: { 
        src: "./src/assets/logo.svg", 
        replacesTitle: true, 
      }, 
      favicon: "./src/assets/logo.svg", 
      social: { 
        github: "https://github.com/Shoaib-Programmer/PyChemicals", 
      }, 
      sidebar: [ 
        { 
          label: "Guides", 
          items: [{ label: "Example Guide", link: "/guides/example/" }], 
        }, 
        { 
          label: "Reference", 
          autogenerate: { directory: "reference" }, 
        }, 
      ], 
      customCss: ["./src/styles/globals.css"], 
      // Removed multi-language options: defaultLocale and locales 
    }), 
  ], 
    markdown: {
    remarkPlugins: [remarkMath],
    rehypePlugins: [rehypeMathjax],
  },
}); 
