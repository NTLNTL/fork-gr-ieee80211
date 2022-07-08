/*
 *
 *     GNU Radio IEEE 802.11a/g/n/ac 2x2
 *     Encoder of 802.11a/g/n/ac 1x1 and 2x2 payload part
 *     Copyright (C) June 1, 2022  Zelin Yun
 *
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Affero General Public License as published
 *     by the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Affero General Public License for more details.
 *
 *     You should have received a copy of the GNU Affero General Public License
 *     along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <gnuradio/io_signature.h>
#include "encode_impl.h"

namespace gr {
  namespace ieee80211 {

    encode::sptr
    encode::make(const std::string& tsb_tag_key)
    {
      return gnuradio::make_block_sptr<encode_impl>(tsb_tag_key
        );
    }


    /*
     * The private constructor
     */
    encode_impl::encode_impl(const std::string& tsb_tag_key)
      : gr::tagged_stream_block("encode",
              gr::io_signature::make(0, 0, 0),
              gr::io_signature::make(2, 2, sizeof(uint8_t)), tsb_tag_key)
    {
      //message_port_register_in(pdu::pdu_port_id());
      d_sEncode = ENCODE_S_IDLE;

      message_port_register_in(pmt::mp("pdus"));
      set_msg_handler(pmt::mp("pdus"), boost::bind(&encode_impl::msgRead, this, _1));
    }

    /*
     * Our virtual destructor.
     */
    encode_impl::~encode_impl()
    {
    }

    void
    encode_impl::msgRead(pmt::pmt_t msg)
    {
      std::cout<<"ieee80211 encode, new msg";
      pmt::pmt_t vector = pmt::cdr(msg);
      int tmpMsgLen = pmt::blob_length(vector);
      size_t tmpOffset(0);
      const uint8_t *tmpPkt = (const uint8_t *)pmt::uniform_vector_elements(vector, tmpOffset);
      uint16_t tmpLen;
      // 1B format, 1B mcs, 1B nss, 2B len, total 5B, len is 0 then NDP
      if((tmpMsgLen < 5) || (tmpMsgLen > DECODE_D_MAX)){
        return;
      }
      //memcpy(d_msg, tmpPkt, tmpMsgLen);
      tmpLen = (((uint16_t)tmpPkt[3])<<8  + (uint16_t)tmpPkt[4]);
      formatToModSu(&d_m, (int)tmpPkt[0], (int)tmpPkt[1], (int)tmpPkt[2], (int)tmpLen);
      if(d_m.format == C8P_F_L)
      {
        // legacy
        legacySigBitsGen(d_legacySig, d_legacySigCoded, d_m.mcs, d_m.len);
        procIntelLegacyBpsk(d_legacySigCoded, d_legacySigInted);

        uint8_t* tmpDataP = d_dataBits;
        memset(tmpDataP, 0, 16);
        tmpDataP += 16;
        for(int i=0;i<d_m.len;i++)
        {
          for(int j=0;j<8;j++)
          {
            tmpDataP[j] = (tmpPkt[i] >> j);
          }
          tmpDataP += 8;
        }
        // tail
        memset(tmpDataP, 0, 6);
        tmpDataP += 6;
        // pad
        memset(tmpDataP, 0, (d_m.nSym * d_m.nDBPS - 22 - d_m.len*8));
      }
      else if(d_m.format == C8P_F_VHT)
      {
        // vht
        vhtSigABitsGenSU(d_vhtSigA, d_vhtSigACoded, &d_m);
        procIntelLegacyBpsk(&d_vhtSigACoded[0], &d_vhtSigAInted[0]);
        procIntelLegacyBpsk(&d_vhtSigACoded[48], &d_vhtSigAInted[48]);
        vhtSigB20BitsGenSU(d_vhtSigB20, d_vhtSigB20Coded, d_vhtSigBCrc8, &d_m);
        procIntelVhtB20(d_vhtSigB20Coded, d_vhtSigB20Inted);

        int tmpPsduLen = (d_m.nSym * d_m.nDBPS - 16 - 6)/8;
        // legacy training 16, legacy sig 4, vhtsiga 8, vht training 4+4n, vhtsigb, payload, no short GI
        int tmpTxTime = 20 + 8 + 4 + d_m.nLTF * 4 + 4 + d_m.nSym * 4;
        int tmpLegacyLen = ((tmpTxTime - 20) / 4 + (((tmpTxTime - 20) % 4) != 0)) * 3 - 3;
        legacySigBitsGen(d_legacySig, d_legacySigCoded, 0, tmpLegacyLen);
        procIntelLegacyBpsk(d_legacySigCoded, d_legacySigInted);

        uint8_t* tmpDataP = d_dataBits;
        // 7 scrambler init, 1 reserved
        memset(tmpDataP, 0, 8);
        tmpDataP += 8;
        // 8 sig b crc8
        memcpy(tmpDataP, d_vhtSigBCrc8, 8);
        tmpDataP += 8;
        // data
        for(int i=0;i<d_m.len;i++)
        {
          for(int j=0;j<8;j++)
          {
            tmpDataP[j] = (tmpPkt[i] >> j);
          }
          tmpDataP += 8;
        }
        // EOF subframe padding
        for(int i=0;i<(tmpPsduLen/4);i++)
        {
          memcpy(tmpDataP, EOF_PAD_SUBFRAME, 32);
          tmpDataP += 32;
        }
        // EOF octect padding
        for(int i=0;i<(tmpPsduLen%4);i++)
        {
          memset(tmpDataP, 0, 8);
          tmpDataP += 8;
        }
        // tail pading, all 0, includes tail bits, when scrambling, do not scramble tail
        memset(tmpDataP, 0, (d_m.nSym * d_m.nDBPS - tmpPsduLen*8 - 16));
      }
      else
      {
        // ht
        htSigBitsGen(d_htSig, d_htSigCoded, &d_m);
        procIntelLegacyBpsk(&d_htSigCoded[0], &d_htSigInted[0]);
        procIntelLegacyBpsk(&d_htSigCoded[48], &d_htSigInted[48]);
        // legacy training and sig 20, htsig 8, ht training 4+4n, payload, no short GI
        int tmpTxTime = 20 + 8 + 4 + d_m.nLTF * 4 + d_m.nSym * 4;
        int tmpLegacyLen = ((tmpTxTime - 20) / 4 + (((tmpTxTime - 20) % 4) != 0)) * 3 - 3;
        legacySigBitsGen(d_legacySig, d_legacySigCoded, 0, tmpLegacyLen);
        procIntelLegacyBpsk(d_legacySigCoded, d_legacySigInted);

        uint8_t* tmpDataP = d_dataBits;
        // service
        memset(tmpDataP, 0, 16);
        tmpDataP += 16;
        // data
        for(int i=0;i<d_m.len;i++)
        {
          for(int j=0;j<8;j++)
          {
            tmpDataP[j] = (tmpPkt[i] >> j);
          }
          tmpDataP += 8;
        }
        // tail
        memset(tmpDataP, 0, 6);
        tmpDataP += 6;
        // pad
        memset(tmpDataP, 0, (d_m.nSym * d_m.nDBPS - 22 - d_m.len*8));
      }
      d_sEncode = ENCODE_S_SCEDULE;
    }

    int
    encode_impl::calculate_output_stream_length(const gr_vector_int &ninput_items)
    {
      if(d_sEncode == ENCODE_S_SCEDULE)
      {
        std::cout<<"schedule in calculate"<<std::endl;
        d_nChipsGen = d_m.nSym * d_m.nSD;      // gen payload part qam chips
        d_nChipsGenProcd = 0;
        d_sEncode = ENCODE_S_ENCODE;
      }
      return d_nChipsGen;
    }

    int
    encode_impl::work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      //auto in = static_cast<const input_type*>(input_items[0]);
      uint8_t* outChips1 = static_cast<uint8_t*>(output_items[0]);
      uint8_t* outChips2 = static_cast<uint8_t*>(output_items[1]);
      d_nGen = noutput_items;

      switch(d_sEncode)
      {
        case ENCODE_S_IDLE:
        {
          std::cout<<"idle"<<std::endl;
          return 0;
        }

        case ENCODE_S_SCEDULE:
        {
          std::cout<<"schedule in work"<<std::endl;
          return 0;
        }

        case ENCODE_S_ENCODE:
        {
          std::cout<<"encode and gen tag"<<std::endl;
          // scrambling
          if(d_m.format == C8P_F_VHT)
          {
            scramEncoder(d_dataBits, d_scramBits, (d_m.nSym * d_m.nDBPS - 6), 97);
            memset(&d_scramBits[d_m.nSym * d_m.nDBPS - 6], 0, 6);
          }
          else
          {
            scramEncoder(d_dataBits, d_scramBits, (d_m.nSym * d_m.nDBPS), 97);
            memset(&d_scramBits[d_m.len * 8 + 16], 0, 6);
          }
          // binary convolutional coding
          bccEncoder(d_scramBits, d_convlBits, d_m.nSym * d_m.nDBPS);
          // puncturing
          punctEncoder(d_scramBits, d_punctBits, d_m.nSym * d_m.nDBPS * 2, &d_m);
          // interleave and convert to qam chips
          if(d_m.nSS == 1)
          {
            if(d_m.format == C8P_F_L)
            {
              for(int i=0;i<d_m.nSym;i++)
              {
                procInterLegacy(&d_punctBits[i*d_m.nCBPS], &d_IntedBits1[i*d_m.nCBPS], &d_m);
              }
            }
            else
            {
              for(int i=0;i<d_m.nSym;i++)
              {
                procInterNonLegacy(&d_punctBits[i*d_m.nCBPS], &d_IntedBits1[i*d_m.nCBPS], &d_m);
              }
            }
            bitsToChips(d_IntedBits1, d_qamChips1, d_m.nSym * d_m.nCBPS, &d_m);
          }
          else
          {
            // stream parser first
            streamParser2(d_punctBits, d_parsdBits1, d_parsdBits2, d_m.nSym * d_m.nCBPS, &d_m);
            for(int i=0;i<d_m.nSym;i++)
            {
              procInterNonLegacy(&d_parsdBits1[i*d_m.nCBPSS], &d_IntedBits1[i*d_m.nCBPSS], &d_m);
              procInterNonLegacy(&d_parsdBits2[i*d_m.nCBPSS], &d_IntedBits2[i*d_m.nCBPSS], &d_m);
            }
            bitsToChips(d_IntedBits1, d_qamChips1, d_m.nSym * d_m.nCBPSS, &d_m);
            bitsToChips(d_IntedBits2, d_qamChips2, d_m.nSym * d_m.nCBPSS, &d_m);
          }

          // gen tag
          d_tagLegacyBits.clear();
          d_tagLegacyBits.reserve(48);
          for(int i=0;i<48;i++)
          {
            d_tagLegacyBits.push_back(d_legacySigInted[i]);
          }
          pmt::pmt_t dict = pmt::make_dict();
          dict = pmt::dict_add(dict, pmt::mp("format"), pmt::from_long(d_m.format));
          dict = pmt::dict_add(dict, pmt::mp("nss"), pmt::from_long(d_m.format));
          dict = pmt::dict_add(dict, pmt::mp("nsym"), pmt::from_long(d_m.nSym));
          dict = pmt::dict_add(dict, pmt::mp("lsig"), pmt::init_u8vector(d_tagLegacyBits.size(), d_tagLegacyBits));
          if(d_m.format == C8P_F_HT)
          {
            d_tagHtBits.clear();
            d_tagHtBits.reserve(96);
            for(int i=0;i<96;i++)
            {
              d_tagHtBits.push_back(d_htSigInted[i]);
            }
            dict = pmt::dict_add(dict, pmt::mp("htsig"), pmt::init_u8vector(d_tagHtBits.size(), d_tagHtBits));
          }
          else if(d_m.format == C8P_F_VHT)
          {
            d_tagVhtABits.clear();
            d_tagVhtABits.reserve(96);
            d_tagVhtB20Bits.clear();
            d_tagVhtB20Bits.reserve(52);
            for(int i=0;i<96;i++)
            {
              d_tagVhtABits.push_back(d_vhtSigAInted[i]);
            }
            for(int i=0;i<52;i++)
            {
              d_tagVhtB20Bits.push_back(d_vhtSigB20Inted[i]);
            }
            dict = pmt::dict_add(dict, pmt::mp("vhtsiga"), pmt::init_u8vector(d_tagVhtABits.size(), d_tagVhtABits));
            dict = pmt::dict_add(dict, pmt::mp("vhtsigb"), pmt::init_u8vector(d_tagVhtB20Bits.size(), d_tagVhtB20Bits));
          }
          pmt::pmt_t pairs = pmt::dict_items(dict);
          for (int i = 0; i < pmt::length(pairs); i++) {
              pmt::pmt_t pair = pmt::nth(i, pairs);
              add_item_tag(0,                   // output port index
                            nitems_written(0),  // output sample index
                            pmt::car(pair),     
                            pmt::cdr(pair),
                            alias_pmt());
          }

          d_sEncode = ENCODE_S_COPY;
          return 0;
        }

        case ENCODE_S_COPY:
        {
          std::cout<<"copy"<<std::endl;
          int o1 = 0;
          while((o1 + d_m.nSD) < d_nGen)
          {
            if(d_m.nSS == 1)
            {
              memcpy(&outChips1[o1], &d_qamChips1[d_nChipsGenProcd], d_m.nSD);
              memset(&outChips2[o1], 0, d_m.nSD);
            }
            else
            {
              memcpy(&outChips1[o1], &d_qamChips1[d_nChipsGenProcd], d_m.nSD);
              memcpy(&outChips2[o1], &d_qamChips2[d_nChipsGenProcd], d_m.nSD);
            }
            o1 += d_m.nSD;
            d_nChipsGenProcd += d_m.nSD;
            if(d_nChipsGenProcd >= d_nChipsGen)
            {
              std::cout<<"copy done"<<std::endl;
              d_sEncode = ENCODE_S_IDLE;
              break;
            }
          }
          return o1;
        }
      }

      // Tell runtime system how many output items we produced.
      d_sEncode = ENCODE_S_IDLE;
      return 0;
    }

  } /* namespace ieee80211 */
} /* namespace gr */
