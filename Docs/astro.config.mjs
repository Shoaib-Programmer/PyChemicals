// @ts-check 
import { defineConfig } from "astro/config"; 
import starlight from "@astrojs/starlight"; 
import remarkMath from "remark-math";
import rehypeMathjax from "rehype-mathjax";
 
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
      favicon: "/favicon.ico", 
      social: { 
        github: "https://github.com/Shoaib-Programmer/PyChemicals", 
      }, 
      sidebar: [ 
        { 
          label: "Guides", 
          items: [{ label: "Getting Started", link: "guides/quickstart" }], 
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
