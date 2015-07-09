/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2014 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** GNU GENERAL PUBLIC LICENSE Version 3
 **
 ** You should have received a copy of the GNU GENERAL PUBLIC LICENSE Version 3
 ** along with this program. If not, see
 ** <https://github.com/illumina/licenses/>.
 **
 ** \file XmlWriter.cpp
 **
 ** Helper class for working with xml
 **
 ** \author Roman Petrovski
 **/

#include "xml/XmlWriter.hh"


namespace isaac
{
namespace xml
{

int XmlWriter::xmlOutputWriteCallback(void *context, const char *buffer, int len)
{
    std::ostream* pStream = reinterpret_cast<std::ostream*>(context);
    if (!*pStream && !pStream->eof())
    {
//        ISAAC_THREAD_CERR << "stream is not open" << std::endl;
        return -1;
    }

    if (!pStream->write(buffer, len))
    {
//        ISAAC_THREAD_CERR << "stream write failed" << std::endl;
        return -1;
    }
//    ISAAC_THREAD_CERR << "stream written " << len << std::endl;
    return len;
}

XmlWriter::XmlWriter(std::ostream &os) :
    os_(os),
    xmlWriter_(xmlNewTextWriter(xmlOutputBufferCreateIO(xmlOutputWriteCallback, 0, &os_, 0)))
{
    if (!xmlWriter_)
    {
        BOOST_THROW_EXCEPTION(XmlWriterException("xmlNewTextWriter failed"));
    }

    int written = xmlTextWriterStartDocument(xmlWriter_, 0, 0, 0);
    if (-1 == written)
    {
        BOOST_THROW_EXCEPTION(XmlWriterException(std::string("xmlTextWriterStartDocument returned -1")));
    }

    if (-1 == xmlTextWriterSetIndent(xmlWriter_, 1))
    {
        BOOST_THROW_EXCEPTION(XmlWriterException(std::string("xmlTextWriterSetIndent returned -1")));
    }

    if (-1 == xmlTextWriterSetIndentString(xmlWriter_, BAD_CAST "  "))
    {
        BOOST_THROW_EXCEPTION(XmlWriterException(std::string("xmlTextWriterSetIndentString returned -1")));
    }
}

XmlWriter::~XmlWriter()
{
    xmlFreeTextWriter(xmlWriter_);
}

void XmlWriter::close()
{
    int written = xmlTextWriterEndDocument(xmlWriter_);
    if (-1 == written)
    {
          BOOST_THROW_EXCEPTION(XmlWriterException(std::string("xmlTextWriterEndDocument returned -1")));
    }
}

XmlWriter &XmlWriter::startElement(const char *name)
{
    const int written = xmlTextWriterStartElement(xmlWriter_, BAD_CAST name);
    if (-1 == written)
    {
        BOOST_THROW_EXCEPTION(XmlWriterException(std::string("xmlTextWriterStartElement returned -1 for element ") + name));
    }
    return *this;
}

XmlWriter &XmlWriter::endElement()
{
    const int written = xmlTextWriterEndElement(xmlWriter_);
    if (-1 == written)
    {
        BOOST_THROW_EXCEPTION(XmlWriterException(std::string("xmlTextWriterEndElement returned -1")));
    }
    return *this;
}

XmlWriter &XmlWriter::writeText(const char *text)
{
    const int written = xmlTextWriterWriteRaw(xmlWriter_, BAD_CAST text);
    if (-1 == written)
    {
        BOOST_THROW_EXCEPTION(XmlWriterException(std::string("xmlTextWriterWriteRaw returned -1 ")));
    }
    return *this;
}


} // namespace xml
} // namespace isaac
