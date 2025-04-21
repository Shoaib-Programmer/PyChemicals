// @ts-check
import { defineConfig } from "astro/config";
import starlight from "@astrojs/starlight";
import remarkMath from "remark-math";
import rehypeMathjax from "rehype-mathjax";

// Import the logo so Vite processes it.
import logo from "./src/assets/logo.svg";

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
        src: logo,
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
    }),
  ],
  markdown: {
    remarkPlugins: [remarkMath],
    rehypePlugins: [rehypeMathjax],
  },
});