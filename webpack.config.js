let webpack = require("webpack"),
    HtmlWebpackPlugin = require("html-webpack-plugin"),
    CopyWebpackPlugin = require("copy-webpack-plugin"),
    CleanWebpackPlugin = require("clean-webpack-plugin"),
    path = require("path");

let src = path.join(__dirname, "/src/"),
    dist = path.join(__dirname, "/dist/");

module.exports = {
    entry: "./src/main.js",
    output: {
        path: dist,
        filename: 'bundle.js',
    },
    module: {
        loaders: [
            {
                test: /\.js$/,
                exclude: /node_modules/,
                loader: "eslint-loader",
                options: {
                    emitter: true,
                }
            },
            {
                test: /\.css$/,
                use: [
                    {loader: "style-loader"},
                    {loader: "css-loader", options: {modules: true}}
                ],
            }
        ]
    },
    plugins: [
        new HtmlWebpackPlugin({
            hash: true,
            template: "./src/index.html",
            filename: "index.html",
            title: "taipei-metromap"
        }),
        new CopyWebpackPlugin([
            {from: `${src}/geodata`, to: `${dist}/geodata`},
            {from: `${src}/svg`, to: `${dist}/svg`},
        ]),
        new CleanWebpackPlugin([dist], {
            "verbose": true,
            "exclude": [],
        })
    ],
    resolve: {
        extensions: [".js", ".json"]
    }
};